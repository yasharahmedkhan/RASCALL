from flask import Flask, request, jsonify, render_template
from flask_sqlalchemy import SQLAlchemy
import os
import json
import pandas as pd
from rascall import analysis
import subprocess
from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.layouts import column
from bokeh.models import Tabs, Panel, Legend, ColumnDataSource
from bokeh.palettes import Spectral10


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
    return render_template('homepage.html')

#@app.route('/plot')
#def plot_page():
 #   return render_template('plot.html')

@app.route('/search')
def search_page():
    return render_template('search.html')

@app.route('/list')
def list_page():
    return render_template('list.html')

@app.route('/documentation')
def documentation():
    return render_template('documentation.html')

@app.route('/run_rascall', methods=['GET'])
def run_rascall():
    # Get all query parameters
    params = request.args.to_dict()
    
    # Prepare the command
    command = ['./rascall_list']
    for key, value in params.items():
        if value:  # Only add non-empty parameters
            command.extend([f'--{key}', value])
    
    app.logger.info(f"Executing command: {' '.join(command)}")
    
    try:
        # Run the command and capture output
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        
        # Log the output for debugging
        app.logger.info(f"Command output: {result.stdout}")
        
        if result.stdout.strip():
            return jsonify({"output": result.stdout})
        else:
            return jsonify({"output": "No results found."})
    except subprocess.CalledProcessError as e:
        app.logger.error(f"Error executing command: {str(e)}")
        app.logger.error(f"Error output: {e.stderr}")
        return jsonify({"error": f"Error executing command: {e.stderr}"}), 500
    except Exception as e:
        app.logger.error(f"Unexpected error: {str(e)}")
        return jsonify({"error": f"An unexpected error occurred: {str(e)}"}), 500
    
@app.route('/plot', methods=['GET', 'POST'])
def plot_page():
    if request.method == 'POST':
        fg = request.form.get('fg')
        mol = request.form.get('mol')
        mf = request.form.get('mf', 'all')
        
        plot, tab_data = generate_tabbed_plots(fg, mf, mol)
        
        if plot:
            script, div = components(plot)
            return render_template('plot_results.html', 
                                script=script, 
                                div=div, 
                                tab_data=json.dumps(tab_data))
        else:
            return render_template('plot.html', error="No data found for the given input.")
    
    return render_template('plot.html')

@app.route('/plot_functional/<functional_group>')
def plot_functional(functional_group):
    plot, tab_data = generate_tabbed_plots(functional_group, 'all', None)
    if plot:
        script, div = components(plot)
        return render_template('plot_results.html', 
                            script=script, 
                            div=div, 
                            tab_data=json.dumps(tab_data))
    else:
        return render_template('plot.html', error=f"No data found for functional group: {functional_group}")

def generate_tabbed_plots(fg, mf, mol):
    data = analysis.plot(fg, mf, mol, return_data=True)
    
    if not data:
        return None, {}
    
    tabs = []
    tab_data = {}  # Store data for each tab
    
    for molecule, molecule_data in data.items():
        p = figure(title=f"Spectrum for {molecule}", 
                  x_axis_label='Wavenumbers (cm^-1)', 
                  y_axis_label='Intensity',
                  width=800, height=400)
        
        legend_items = []
        colors = Spectral10
        color_index = 0
        
        # Store functional groups and their colors for this molecule
        tab_data[molecule] = {
            'functional_groups': [],
            'colors': [],
            'is_first': len(tab_data) == 0  # Mark first tab
        }
        
        for functional, points in molecule_data['rascall'].items():
            if points:
                color = colors[color_index]
                x, y = zip(*points)
                source = ColumnDataSource(data=dict(x=x, y=y))
                
                stems = p.segment(x0='x', y0=0, x1='x', y1='y', color=color, source=source)
                circles = p.circle('x', 'y', size=5, color=color, source=source)
                
                legend_items.append((functional, [stems, circles]))
                
                # Store functional group data
                tab_data[molecule]['functional_groups'].append(functional)
                tab_data[molecule]['colors'].append(color)
                
                color_index = (color_index + 1) % len(colors)
        
        if 'nist' in molecule_data and molecule_data['nist']:
            x, y = zip(*molecule_data['nist'])
            line = p.line(x, y, line_color='black', line_width=1.5)
            legend_items.append(('NIST', [line]))
        
        legend = Legend(items=legend_items, 
                       location="center_right", 
                       orientation="vertical", 
                       label_text_font_size="8pt")
        legend.click_policy = "hide"
        
        p.add_layout(legend, 'right')
        p.legend.visible = True
        
        tab = Panel(child=p, title=molecule)
        tabs.append(tab)
    
    return Tabs(tabs=tabs), tab_data


# @app.route('/plot_rascall', methods=['GET'])
# def plot_rascall():
#     fg = request.args.get('fg')
#     mol = request.args.get('mol')
#     mf = request.args.get('mf', 'all')
#     aw = request.args.get('aw')
    
#     try:
#         plots = analysis.plot(fg, mf, mol, aw)
        
#         tabs = []
#         for i, (fig, molecule_name) in enumerate(zip(plots, plots.molecule_names)):
#             panel = Panel(child=fig, title=molecule_name)
#             tabs.append(panel)
        
#         tabs = Tabs(tabs=tabs)
        
#         # Wrap the plot into a json_item with an id
#         plot_json = json_item(tabs, "bokeh-plot")
#         return jsonify(plot_json)
#     except Exception as e:
#         return jsonify({"error": str(e)}), 500


@app.route('/search_results', methods=['GET'])
def search():
    query = request.args.get('query')
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
        db.create_all()
        populate_db()
    app.run(debug=True, port=5001)


