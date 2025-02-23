from rdkit import Chem
from rdkit.Chem import Draw
import io
import pandas as pd
import os
from functools import lru_cache
import base64


class MoleculeVisualizer:
    def __init__(self, data_dir='data'):
        self.data_dir = data_dir

    def process_request(self, element_id, filename):
        """process request data"""
        try:
            # generate compound image
            img_io = self.generate_structure_image(element_id, filename)
            # transform image to base64 string
            img_io.seek(0)
            img_base64 = base64.b64encode(img_io.getvalue()).decode()
            return {
                'success': True,
                'data': img_base64,
                'message': 'Image generated successfully'
            }
        except ValueError as e:
            return {
                'success': False,
                'error': str(e)
            }
        except FileNotFoundError as e:
            return {
                'success': False,
                'error': f'File not found: {filename}'
            }
        except Exception as e:
            return {
                'success': False,
                'error': str(e)
            }

    def generate_structure_image(self, element_id, filename):
        """generate structure image of a molecule"""
        try:
            # validate if the file is csv
            if not filename.endswith('.csv'):
                raise ValueError('Only CSV files are supported')
            # search SMILES
            smiles = self._find_smiles_by_id_and_file(element_id, filename)
            if not smiles:
                raise ValueError(f'SMILES not found for ID: {element_id}')
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f'Invalid SMILES: {smiles}')

            #### modify image size

            img = Draw.MolToImage(mol, size=(300, 300), bgcolor='white')
            # convert image to byte stream
            img_io = io.BytesIO()
            img.save(img_io, format='PNG')
            img_io.seek(0)
            return img_io
        except Exception as e:
            raise Exception(f'Error generating structure: {str(e)}')

    @lru_cache(maxsize=1000)
    def _find_smiles_by_id_and_file(self, element_id, filename):
        """find smiles by id and filename"""
        try:
            file_path = os.path.join(self.data_dir, filename)
            df = pd.read_csv(file_path,dtype={'cmpd_id': str})

            cmpd_id_str = f"cmpd_{element_id}"


            result = df[df['cmpd_id'] == cmpd_id_str]
            if result.empty:
                print(f'No matching cmpd_id found for {element_id} in element_id')
                return None
            return result['SMILES'].iloc[0]
        except FileNotFoundError:
            raise FileNotFoundError(f'File not found: {filename}')
        except Exception as e:
            raise Exception(f'Error reading file {filename}: {str(e)}')