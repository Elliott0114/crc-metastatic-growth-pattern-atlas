"""Small offline tests for the public workflow's failure and data contracts."""
from pathlib import Path
import os
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproduce import cells_equal, clean_env, validate_files, verify_provenance
from workflow_io import sha256, read_tsv, read_workbook, write_workbook


class Contracts(unittest.TestCase):
    def test_author_annotation_is_external_only(self):
        root = Path(__file__).resolve().parents[1]
        name = 'visium_liver_meta_after_qc.tsv.gz'
        rows = [r for r in read_tsv(root / 'metadata/analysis_inputs.tsv') if r['path'].endswith(name)]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['location'], 'external')
        self.assertEqual(rows[0]['path'], rows[0]['workspace_path'])
        self.assertIn('authors', rows[0]['source'])
        self.assertFalse((root / 'data/inputs' / rows[0]['workspace_path']).exists())

    def test_missing_and_modified_inputs_stop(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / 'counts.tsv'
            path.write_text('gene\tS1\nA\t2\n', encoding='utf-8')
            row = dict(path=path.name, bytes=path.stat().st_size, sha256=sha256(path), source='source accession')
            self.assertEqual(validate_files([row], root, ''), [])
            path.write_text('gene\tS1\nA\t3\n', encoding='utf-8')
            self.assertIn('Checksum mismatch', validate_files([row], root, '')[0])
            path.unlink()
            self.assertIn('source accession', validate_files([row], root, '')[0])

    def test_unsafe_input_path(self):
        self.assertIn('Unsafe', validate_files([dict(path='../elsewhere')], Path('.'), '')[0])

    def test_discrete_order_missingness_and_float_tolerance(self):
        self.assertTrue(cells_equal([['A', 1.0000001, None, True]], [['A', 1, None, True]]))
        for value in [[['B', 1, None, True]], [['A', 1, 'NA', True]], [['A', 1, None, 1]], [['A', 1.001, None, True]]]:
            self.assertFalse(cells_equal(value, [['A', 1, None, True]]))
        self.assertFalse(cells_equal([['A'], ['B']], [['B'], ['A']]))
        self.assertFalse(cells_equal([[1000001]], [[1000000]]))

    def test_workbook_roundtrip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'book.xlsx'
            sheets = [('Data', [['id', 'n', 'flag', 'missing'], ['A', 2.5, True, None]]), ('Notes', [['NA']])]
            write_workbook(path, sheets)
            self.assertEqual(read_workbook(path), sheets)

    def test_direct_connections(self):
        with patch.dict(os.environ, {'HTTPS_PROXY': 'http://invalid', 'all_proxy': 'http://invalid'}):
            self.assertNotIn('HTTPS_PROXY', clean_env())
            self.assertNotIn('all_proxy', clean_env())

    def test_source_index_cannot_omit_a_declared_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'metadata').mkdir()
            sources = ['a.tsv', 'b.tsv']
            recipes = [dict(table=1, kind='data', sheet='Data', sources=sources)]
            (root / 'metadata/workbook_recipes.json').write_text(json.dumps(recipes), encoding='utf-8')
            rows = [['sheet', 'source', 'sha256']]
            for source in sources:
                path = root / source
                path.write_text('gene\tvalue\nA\t1\n', encoding='utf-8')
                rows.append(['Data', source, sha256(path)])
            recipe = dict(table=1, sheet='Added_source_index')
            with patch('reproduce.ROOT', root):
                self.assertTrue(verify_provenance(recipe, rows, root))
                self.assertFalse(verify_provenance(recipe, rows[:-1], root))


if __name__ == '__main__':
    unittest.main()
