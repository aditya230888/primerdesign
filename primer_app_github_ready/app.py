from flask import Flask, request, render_template, Response
from Bio import Entrez, SeqIO
import primer3
import csv
import io

app = Flask(__name__)
Entrez.email = "your_email@example.com"

PRESET_QUERIES = {
    "Mitochondrial DNA Haplotype": "human mitochondrial HVR1",
    "HVR1": "human mitochondrial HVR1",
    "HVR2": "human mitochondrial HVR2",
    "CYTB": "human mitochondrial CYTB gene",
    "D-loop": "human mitochondrial D-loop",
    "12S rRNA": "human mitochondrial 12S rRNA"
}

def score_primer(primer_seq):
    if not primer_seq:
        return {'score': 0, 'tm': 0, 'gc': 0, 'warnings': ['No sequence']}
    tm = primer3.calcTm(primer_seq)
    gc = 100.0 * (primer_seq.count('G') + primer_seq.count('C')) / len(primer_seq)
    score = 100
    warnings = []
    if len(primer_seq) < 18 or len(primer_seq) > 25:
        score -= 20
        warnings.append("Unusual length")
    if tm < 55 or tm > 65:
        score -= 20
        warnings.append("Suboptimal Tm")
    if gc < 40 or gc > 60:
        score -= 20
        warnings.append("GC content out of range")
    if primer3.calcHomodimer(primer_seq).dg < -5:
        score -= 20
        warnings.append("Self-dimer risk")
    if primer3.calcHairpin(primer_seq).dg < -5:
        score -= 20
        warnings.append("Hairpin risk")
    return {
        'score': max(score, 0),
        'tm': round(tm, 2),
        'gc': round(gc, 2),
        'warnings': warnings
    }

def fetch_sequence_ncbi(query):
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    record = Entrez.read(handle)
    if record["IdList"]:
        seq_id = record["IdList"][0]
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        return str(seq_record.seq)
    else:
        return None

def design_primers(sequence):
    result = primer3.bindings.designPrimers(
        {
            'SEQUENCE_TEMPLATE': sequence
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PRODUCT_SIZE_RANGE': [[80, 150]],
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_NUM_RETURN': 5
        }
    )
    if result['PRIMER_LEFT_NUM_RETURNED'] == 0:
        raise ValueError("Primer3 could not design primers. Try a different sequence.")
    primers = []
    for i in range(min(result['PRIMER_LEFT_NUM_RETURNED'], 5)):
        fwd = result.get(f'PRIMER_LEFT_{i}_SEQUENCE', '')
        rev = result.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')
        probe = result.get(f'PRIMER_INTERNAL_{i}_SEQUENCE', '')
        primers.append({
            'forward': fwd,
            'reverse': rev,
            'probe': probe,
            'forward_score': score_primer(fwd),
            'reverse_score': score_primer(rev),
            'probe_score': score_primer(probe)
        })
    return primers

@app.route('/', methods=['GET', 'POST'])
def index():
    primers = None
    error = None
    if request.method == 'POST':
        user_input = request.form.get('pathogen') or request.form.get('dropdown')
        query = PRESET_QUERIES.get(user_input.strip(), user_input.strip())
        try:
            sequence = fetch_sequence_ncbi(query)
            if sequence:
                primers = design_primers(sequence)
            else:
                error = "No sequence found for query."
        except Exception as e:
            error = f"Error: {str(e)}"
    return render_template("index.html", primers=primers, error=error)

@app.route('/download', methods=['POST'])
def download_csv():
    query = request.form['query']
    sequence = fetch_sequence_ncbi(PRESET_QUERIES.get(query.strip(), query.strip()))
    primers = design_primers(sequence)
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(['Forward Primer', 'Reverse Primer', 'Probe', 'Tm F', 'GC F', 'Tm R', 'GC R', 'Tm P', 'GC P'])
    for p in primers:
        writer.writerow([
            p['forward'], p['reverse'], p['probe'],
            p['forward_score']['tm'], p['forward_score']['gc'],
            p['reverse_score']['tm'], p['reverse_score']['gc'],
            p['probe_score']['tm'], p['probe_score']['gc']
        ])
    return Response(output.getvalue(), mimetype='text/csv',
                    headers={"Content-disposition": "attachment; filename=primers.csv"})

if __name__ == '__main__':
    app.run(debug=True)
