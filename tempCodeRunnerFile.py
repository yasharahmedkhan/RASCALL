from flask import Flask, request, jsonify, render_template
from flask_sqlalchemy import SQLAlchemy
import subprocess
import os
import json
import pandas as pd

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///RASCALL_Molecule_Identifiers.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class Chemical(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    smiles = db.Column(db.String, unique=True, nullable=False)
    rdkit_smiles = db.Column(db.String, nullable=False)
    molecular_weight = db.Column(db.Float, nullable=False)
    molecular_formula = db.Column(db.String, nullable=False)
    iupac_name = db.Column(db.String, nullable=False)
    episuite_smiles = db.Column(db.String, nullable=False)
    inchi_code = db.Column(db.String, nullable=False)
    inchi_key = db.Column(db.String, nullable=False)

def populate_db():
    if os.path.exists('RASCALL_Molecule_Identifiers.csv'):
        df = pd.read_csv('RASCALL_Molecule_Identifiers.csv')
        for _, row in df.iterrows():
            if not Chemical.query.filter_by(smiles=row['SMILES']).first():
                chemical = Chemical(
                    smiles=row['SMILES'],
                    rdkit_smiles=row['rdkit_SMILES'],
                    molecular_weight=row['Molecular Weight'],
                    molecular_formula=row['Molecular Formula'],
                    iupac_name=row['IUPAC chemical name'],
                    episuite_smiles=row['EPISUITE-compatible SMILES'],
                    inchi_code=row['InChI Code'],
                    inchi_key=row['InChI Key']
                )
                db.session.add(chemical)
        db.session.commit()
    else:
        print("CSV file not found. Please ensure 'RASCALL_Molecule_Identifiers.csv' exists.")

@app.route('/')
def homepage():
    return render_template('index.html')

@app.route('/search_page')
def search_page():
    return render_template('search.html')

@app.route('/run_rascall', methods=['POST'])
def run_rascall():
    mol = request.json.get('mol')
    if not mol:
        return jsonify({"error": "Molecule is not provided"}), 400
    
    try:
        result = subprocess.run(['./rascall_list', '--mol', mol], capture_output=True, text=True)
        return jsonify({"output": result.stdout})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/plot_rascall', methods=['POST'])
def plot_rascall():
    mol = request.json.get('mol')
    fg = request.json.get('fg')
    mf = request.json.get('mf')
    fw = request.json.get('fw')
    ap = request.json.get('ap')
    if not mol and not fg:
        return jsonify({"error": "No molecule or functional group is provided"}), 400
    
    try:
        command = ['./rascall_plot']
        if mol:
            command.extend(['--mol', mol])
        if fg:
            command.extend(['--fg', fg])
        if mf:
            command.extend(['--mf', mf])
        if fw:
            command.extend(['--fw', fw])
        if ap:
            command.extend(['--ap', ap])

        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            return jsonify({"error": result.stderr}), 500
        
        plot_data = json.loads(result.stdout)
        if 'error' in plot_data:
            return jsonify({"error": plot_data['error']}), 400
        return jsonify(plot_data)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/search', methods=['GET'])
def search():
    query = request.args.get('query')
    print(f"Search query received: {query}")  # Debug statement
    if not query:
        return jsonify({'error': 'No query provided'}), 400

    try:
        query_float = float(query)
        float_condition = (Chemical.molecular_weight == query_float)
    except ValueError:
        float_condition = False

    results = Chemical.query.filter(
        (Chemical.smiles == query) |
        (Chemical.rdkit_smiles == query) |
        (Chemical.molecular_formula == query) |
        (Chemical.iupac_name == query) |
        (Chemical.episuite_smiles == query) |
        (Chemical.inchi_code == query) |
        (Chemical.inchi_key == query) |
        (float_condition if float_condition else False)
    ).all()

    results_data = [{
        'SMILES': chemical.smiles,
        'rdkit_SMILES': chemical.rdkit_smiles,
        'Molecular Weight': chemical.molecular_weight,
        'Molecular Formula': chemical.molecular_formula,
        'IUPAC Name': chemical.iupac_name,
        'EPISUITE-Compatible SMILES': chemical.episuite_smiles,
        'InChI Code': chemical.inchi_code,
        'InChI Key': chemical.inchi_key
    } for chemical in results]

    return jsonify(results_data)

if __name__ == '__main__':
    with app.app_context():
        db.drop_all()
        db.create_all()
        populate_db()

    app.run(debug=True)