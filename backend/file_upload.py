# file_upload.py
from flask import request, jsonify
import os
from molecule_annotate import service  # Import the service instance

# Configure the directory for storing uploaded files
UPLOAD_FOLDER = 'data/'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

def upload_file():
    """Handle file upload"""
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']

    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file:
        try:
            # Save the file
            file_path = os.path.join(UPLOAD_FOLDER, file.filename)
            file.save(file_path)

            # Load compound data
            if service.load_compounds(file_path) is not None and not service.compounds_df.empty:
                return jsonify({
                    'message': 'File uploaded and compounds loaded successfully',
                    'filename': file.filename
                }), 200
            else:
                return jsonify({'error': 'Error loading compounds from file'}), 500
        except Exception as e:
            return jsonify({'error': f'Error processing file: {str(e)}'}), 500