
#!/usr/bin/env python3
from dotenv import load_dotenv
from flask import Flask, jsonify, send_from_directory, send_file, Response
from flask_httpauth import HTTPBasicAuth
from werkzeug.security import generate_password_hash, check_password_hash
import os
import json
import requests

load_dotenv()
app = Flask(__name__)
auth = HTTPBasicAuth()

users = {
    "admin": generate_password_hash(os.getenv('FLASK_PASSWORD'))
}

@auth.verify_password
def verify_password(username, password):
    if username in users and check_password_hash(users.get(username), password):
        return username

@app.before_request
def require_login():
    if not auth.current_user():
        return auth.login_required(lambda: None)()

@app.route('/')
def index():
    return send_file('variant_viewer.html')

@app.route('/<path:filename>')
def static_files(filename):
    return send_from_directory('.', filename)

@app.route('/api/files')
def file_list():
    data_dir = 'data'
    files = []
    
    if os.path.exists(data_dir):
        for filename in os.listdir(data_dir):
            if filename.endswith('.json'):
                filepath = os.path.join(data_dir, filename)
                try:
                    stat = os.stat(filepath)
                    files.append({
                        'name': filename,
                        'size': stat.st_size,
                        'modified': stat.st_mtime * 1000
                    })
                except Exception as e:
                    print(f"Error reading {filename}: {e}")
    
    return jsonify(files)

@app.route('/build/<filename>')
def serve_build(filename):
    return send_from_directory('build', filename)

@app.route('/data/<filename>')
def serve_data(filename):
    return send_from_directory('data', filename)

@app.route('/logos/<filename>')
def serve_logo(filename):
    return send_from_directory('logos', filename)

@app.route('/genomes/<path:filename>')
def serve_genome(filename):
    return send_from_directory('igv-genomes', filename)

@app.route('/gos-dist/<path:filename>')
def serve_gosling(filename):
    return send_from_directory('gos-dist/gosling.js/dist', filename)

@app.route('/api/genomes')
def genome_list():
    genome_file = os.path.join('igv-genomes', 'genomes.json')
    if not os.path.exists(genome_file):
        return jsonify({"error": "genomes.json not found"}), 404
    with open(genome_file, 'r') as f:
        genome_config = json.load(f)

    # Add annotation track if GTF/GFF present
    for genome_id, genome in genome_config.items():
        genome_dir = os.path.join('igv-genomes', genome_id)
        if os.path.isdir(genome_dir):
            # Look for .gtf or .gff file
            for fname in os.listdir(genome_dir):
                if fname.endswith('.gtf') or fname.endswith('.gff'):
                    genome['annotationURL'] = f"/genomes/{genome_id}/{fname}"
                    genome['annotationFormat'] = 'gtf' if fname.endswith('.gtf') else 'gff'
                    break  # Only pick the first found

    return jsonify(genome_config)

@app.route('/api/igv-proxy/<filename>')
def igv_proxy(filename):
    allowed_files = ['igv.min.js', 'igv.css', 'igv.esm.min.js']
    
    if filename not in allowed_files:
        return "File not found", 404
    
    # Check if file exists locally
    local_path = os.path.join('igv-dist', filename)
    if os.path.exists(local_path):
        return send_from_directory('igv-dist', filename)
    
    # Download from CDN if not cached
    cdn_url = f'https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/{filename}'
    
    try:
        response = requests.get(cdn_url)
        if response.status_code == 200:
            # Create igv-dist directory if it doesn't exist
            os.makedirs('igv-dist', exist_ok=True)
            
            # Save to local cache
            with open(local_path, 'wb') as f:
                f.write(response.content)
            
            # Return file with appropriate content type
            if filename.endswith('.js'):
                content_type = 'application/javascript'
            elif filename.endswith('.css'):
                content_type = 'text/css'
            else:
                content_type = 'application/octet-stream'
                
            return Response(response.content, mimetype=content_type)
        else:
            return "Download failed", 500
    except Exception as e:
        print(f"Error downloading {filename}: {e}")
        return "Proxy error", 500

if __name__ == '__main__':
    os.makedirs('build', exist_ok=True)
    os.makedirs('data', exist_ok=True)
    os.makedirs('igv-dist', exist_ok=True)
    os.makedirs('igv-genomes/danRer7', exist_ok=True)
    app.run(host='0.0.0.0', port=5010, debug=True, ssl_context="adhoc")