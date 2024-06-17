# from flask import Flask, request, jsonify, send_from_directory
# import subprocess
# import os
# import json

# app = Flask(__name__)


# @app.route('/')
# def homepage():
#     return send_from_directory(os.path.join(app.root_path, 'static'), 'index.html')

# @app.route('/run_rascall', methods=['POST']) 

# def run_rascall():
#     mol = request.json.get('mol')
#     if not mol:
#         return jsonify({"error": "Molecule is not provided"}), 400
    
#     try:
#         result = subprocess.run(['./rascall_list', '--mol', mol], capture_output=True, text=True)
#         return jsonify({"output": result.stdout})
#     except Exception as e:
#         return jsonify({"error": str(e)}), 500

# @app.route('/plot_rascall', methods=['POST'])
# def plot_rascall():
#     mol = request.json.get('mol')
#     fg = request.json.get('fg')
#     mf = request.json.get('mf')
#     fw = request.json.get('fw')
#     ap = request.json.get('ap')
#     if not mol and not fg:
#         return jsonify({"error": "No moleculte or functional group is provided"}), 400
    
#     try:
#         command = ['./rascall_plot']

#         if mol:
#             command.extend(['--mol', mol])
#         if fg:
#             command.extend(['--fg', fg])
#         if mf:
#             command.extend(['--mf', mf])
#         if fw:
#             command.extend(['--fw', fw])
#         if ap:
#             command.extend(['--ap', ap])

#         result = subprocess.run(command, capture_output=True, text=True)
#         if result.returncode != 0:
#             return jsonify({"error": result.stderr}), 500
        
#         plot_json = result.stdout
#         plot_data = json.loads(plot_json)
#         return jsonify(plot_data)
#     except Exception as e:
#         return jsonify({"error": str(e)}), 500
    
# if __name__ == '__main__':
#     app.run(debug=True)