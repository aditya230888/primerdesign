# qPCR Primer/Probe Designer (Flask App)

This web app allows users to design qPCR primers and hydrolysis probes for genes or pathogens by entering either specific gene names or selecting from known mitochondrial targets. It includes:

- Sequence fetching from NCBI
- Primer/probe design using Primer3
- Quality scoring (Tm, GC%, dimers, hairpins)
- CSV export and clipboard copy
- Support for common mitochondrial targets via dropdown

## 🛠 Installation

```bash
git clone https://github.com/your-username/primer-designer.git
cd primer-designer
pip install -r requirements.txt
```

## 🚀 Run the App

```bash
python app.py
```

Then open your browser at: [http://127.0.0.1:5000](http://127.0.0.1:5000)

## 🌐 Deploy on Render.com

1. Push this code to a GitHub repo
2. Go to [https://render.com](https://render.com)
3. Create a new Web Service from your repo
4. Use the following:
   - **Build Command**: `pip install -r requirements.txt`
   - **Start Command**: `python app.py`
   - **Environment**: Python 3.10+ preferred

## ✅ Example Inputs

- `Chlamydia trachomatis`
- `Plasmodium falciparum 18S rRNA`
- `Mitochondrial DNA Haplotype` → auto-mapped to HVR1
- Mito dropdown: `HVR1`, `HVR2`, `CYTB`, `D-loop`, `12S rRNA`

---

Built by RetroBiotech for internal molecular diagnostic development.
