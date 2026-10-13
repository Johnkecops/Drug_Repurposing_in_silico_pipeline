"""
Shared SMILES reader. Gives every compound one stable ID across pipeline stages.

Stages 2 (library prep), 3 (docking) and 7 (ADMET) must agree on compound IDs,
or stage 8 joins the wrong rows. The ID is the 0-based position among SMILES
records, after skipping blank lines and '#' comments.
"""


def read_smiles(path):
    """Return [(compound_id, smiles), ...] for each record in a .smi file."""
    records = []
    with open(path) as f:
        for line in f:
            text = line.strip()
            if not text or text.startswith('#'):
                continue
            records.append((len(records), text.split()[0]))
    return records
