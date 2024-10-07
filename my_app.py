# from flask import Flask, request, jsonify, render_template
# from flask_sqlalchemy import SQLAlchemy
# import pandas as pd
# import os

# app = Flask(__name__)
# app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///RASCALL_Molecule_Identifiers.db'
# app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
# db = SQLAlchemy(app)

# class Chemical(db.Model):
#     id = db.Column(db.Integer, primary_key=True)
#     smiles = db.Column(db.String, unique=True, nullable=False)
#     rdkit_smiles = db.Column(db.String, nullable=False)
#     molecular_weight = db.Column(db.Float, nullable=False)
#     molecular_formula = db.Column(db.String, nullable=False)
#     iupac_name = db.Column(db.String, nullable=False)
#     episuite_smiles = db.Column(db.String, nullable=False)
#     inchi_code = db.Column(db.String, nullable=False)
#     inchi_key = db.Column(db.String, nullable=False)

# def populate_db():
#     # Ensure the CSV file exists before attempting to read
#     if os.path.exists('RASCALL_Molecule_Identifiers.csv'):
#         df = pd.read_csv('RASCALL_Molecule_Identifiers.csv')
#         for _, row in df.iterrows():
#             if not Chemical.query.filter_by(smiles=row['SMILES']).first():
#                 chemical = Chemical(
#                     smiles=row['SMILES'],
#                     rdkit_smiles=row['rdkit_SMILES'],
#                     molecular_weight=row['Molecular Weight'],
#                     molecular_formula=row['Molecular Formula'],
#                     iupac_name=row['IUPAC chemical name'],
#                     episuite_smiles=row['EPISUITE-compatible SMILES'],
#                     inchi_code=row['InChI Code'],
#                     inchi_key=row['InChI Key']
#                 )
#                 db.session.add(chemical)
#         db.session.commit()
#     else:
#         print("CSV file not found. Please ensure 'RASCALL_Molecule_Identifiers.csv' exists.")

# @app.route('/')
# def index():
#     return render_template('search.html')

# @app.route('/search')
# def search():
#     query = request.args.get('query')
#     if not query:
#         return jsonify({'error': 'No query provided'}), 400

#     # Search by exact match
#     results = Chemical.query.filter(
#         (Chemical.smiles == query) |
#         (Chemical.rdkit_smiles == query) |
#         (Chemical.molecular_formula == query) |
#         (Chemical.iupac_name == query) |
#         (Chemical.episuite_smiles == query) |
#         (Chemical.inchi_code == query) |
#         (Chemical.inchi_key == query) |
#         (Chemical.molecular_weight == float(query) if query.replace('.', '', 1).isdigit() else False)
#     ).all()

#     # Return results in JSON format
#     return jsonify([{
#         'SMILES': chemical.smiles,
#         'rdkit_SMILES': chemical.rdkit_smiles,
#         'Molecular Weight': chemical.molecular_weight,
#         'Molecular Formula': chemical.molecular_formula,
#         'IUPAC Name': chemical.iupac_name,
#         'EPISUITE-Compatible SMILES': chemical.episuite_smiles,
#         'InChI Code': chemical.inchi_code,
#         'InChI Key': chemical.inchi_key
#     } for chemical in results])

# if __name__ == '__main__':
#     # Ensure database is recreated with the correct schema
#     with app.app_context():
#         db.drop_all()
#         db.create_all()
#         populate_db()

#     # Run the Flask app
#     app.run(debug=True)






