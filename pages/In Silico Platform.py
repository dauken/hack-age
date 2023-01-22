import streamlit as st
from stmol import showmol
import py3Dmol
import random
import pandas as pd
import requests
from Bio.PDB import PDBList, PDBParser
import biotite.structure.io as bsio
from PIL import Image

from rdkit import Chem
from rdkit.Chem import Draw


@st.cache(allow_output_mutation=True)
def download_dataset():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv(
        "https://raw.githubusercontent.com/aspuru-guzik-group/chemical_vae/master/models/zinc_properties/250k_rndm_zinc_drugs_clean_3.csv"
    ).dropna()
    return df[['smiles']].sample(50)

url = 'https://passer.smu.edu/api'

def pocket_detection(protein):
    data = {"pdb": '5dkk', "chain": "A"}
    results = requests.post(url, data=data)
    print("5")
    print(results.json())
    print("6")
    pocket_residues = results.json()["1"]["residues"].split(" ")[4:]
    pocket_residues = [eval(i) for i in pocket_residues]
    return pocket_residues

def get_sequence(protein, txt_chain):

    pdbl = PDBList()

    # Download the PDB file from the internet
    pdb_file = pdbl.retrieve_pdb_file(protein, file_format='pdb')

    # Create a PDB parser object
    parser = PDBParser()

    # Load the PDB file
    structure = parser.get_structure(protein, pdb_file)

    # Extract the amino acid sequences from the PDB file
    sequence = ""
    for model in structure:
        for chain in model:
            if chain.get_id() == txt_chain:
                for residue in chain:
                    # Only consider amino acid residues
                    if residue.get_resname() in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
                                                 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]:
                        sequence += residue.get_resname()
    if protein == "7WT1":
        sequence = "MAEPQPPSGGLTDEAALSCCSDADPSTKDFLLQQTMLRVKDPKKSLDFYTRVLGMTLIQKCDFPIMKFSLYFLAYEDKNDIPKEKDEKIAWALSRKATLELTHNWGTEDDETQSYHNGNSDPRGFGHIGIAVPDVYSACKRFEELGVKFVKKPDDGKMKGLAFIQDPDGYWIEILNPNKMATLM"
    st.markdown(sequence)
    return sequence


def get_mutations(sequence, mutation_rate):
    st.markdown("""
    Mutated amino acid sequences for selected protein
    """)
    dna_list = list(sequence)
    for i in range(len(sequence)):
        r = random.random()
        if r < mutation_rate:
            mutation_site = random.randint(0, len(dna_list) - 1)
            print(mutation_site)
            dna_list[mutation_site] = random.choice(list('ARNDCQEGHILKMFPSTWYV'))
            mutated_sequence = ''.join(dna_list)

    st.markdown(mutated_sequence)
    return mutated_sequence

# stmol
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb,'pdb')
    pdbview.setStyle({'cartoon':{'color':'spectrum'}})
    pdbview.setBackgroundColor('white')#('0xeeeeee')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height = 500,width=800)

def protein_folding(sequence):
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    name = sequence[:3] + sequence[-3:]
    pdb_string = response.content.decode('utf-8')

    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)

    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 4)

    # Display protein structure
    st.subheader('Mutated protein structure predicted by ESMFold')
    render_mol(pdb_string)

    # plDDT value is stored in the B-factor field
    st.subheader('plDDT')
    st.write('plDDT is a per-residue estimate of the confidence in prediction on a scale from 0-100.')
    st.info(f'plDDT: {b_value}')

    st.download_button(
        label="Download PDB",
        data=pdb_string,
        file_name='predicted.pdb',
        mime='text/plain',
    )


def main():
    st.sidebar.title('In Silico Directed Evolution Platform')
    prot_str='1BML,1D5M,1D5X,1D5Z,1D6E,1DEE,1E9F,1FC2,1FCC,1G4U,1GZS,1HE1,1HEZ,1HQR,1HXY,1IBX,1JBU,1JWM,1JWS'
    prot_list=prot_str.split(',')
    bcolor = "white"
    # Protein sequence input
    DEFAULT_PDB = "7WT1"
    protein = st.sidebar.text_area('Input sequence', DEFAULT_PDB)
    # Protein sequence input
    DEFAULT_CHAIN = "A"
    txt_chain = st.sidebar.text_area('Input chain', DEFAULT_CHAIN)
    # Protein sequence input
    DEFAULT_SMILES = "CC(=O)C(O)SCC(C(=O)NCC(=O)O)NC(=O)CCC(C(=O)O)N"
    txt_smiles = st.sidebar.text_area('Input SMILES', DEFAULT_SMILES, height=100)

    style = st.sidebar.selectbox('style', ['cartoon', 'line', 'cross', 'stick', 'sphere', 'clicksphere'])
    # residues = pocket_detection(protein)
    xyzview = py3Dmol.view(query='pdb:'+protein)
    xyzview.setStyle({style:{'color':'spectrum'}})
    # xyzview.addStyle({'within': {'distance': 3, 'sel': {'chain': txt_chain, 'resi': residues}}}, {'sphere': {'color': 'red'}})
    # xyzview.addStyle({'within': {'distance': 3, 'sel': {'chain': txt_chain}}}, {'sphere': {'color': 'red'}})
    xyzview.setBackgroundColor(bcolor)
    xyzview.spin(True)

    if "button_clicked1" not in st.session_state:
        st.session_state.button_clicked1 = False
    if "button_clicked2" not in st.session_state:
        st.session_state.button_clicked2 = False
    if "button_clicked3" not in st.session_state:
        st.session_state.button_clicked3 = False

    def callback1():
        st.session_state.button_clicked1 = True
    def callback2():
        st.session_state.button_clicked2 = True
    def callback3():
        st.session_state.button_clicked3 = True

    button1 = st.sidebar.button('Start project', on_click = callback1())

    if (button1 or st.session_state.button_clicked1):
        st.subheader('Selected protein visualization')
        showmol(xyzview, height=500, width=800)
        # render_mol(protein)
        st.markdown("""<hr style="height:10px;border:none;color:#333;background-color:#333;" /> """, unsafe_allow_html=True)
        st.subheader('Selected small molecule 2D visualization')
        mol = Chem.MolFromSmiles(txt_smiles)
        img = Draw.MolToImage(mol)
        st.image(img, width=250)

        st.subheader('Amino acid sequence for selected protein')
        sequence = get_sequence(protein, txt_chain)
        st.markdown("""<hr style="height:10px;border:none;color:#333;background-color:#333;" /> """, unsafe_allow_html=True)
        button2 = st.button('Generate mutations', on_click = callback2())
        if (button2 or st.session_state.button_clicked2):
            st.subheader('Amino acid sequence for mutated protein')
            mut8_sequence = get_mutations(sequence, 0.03)
            st.markdown("""<hr style="height:10px;border:none;color:#333;background-color:#333;" /> """,
                        unsafe_allow_html=True)
            protein_folding(mut8_sequence)
            st.markdown("""<hr style="height:10px;border:none;color:#333;background-color:#333;" /> """,
                        unsafe_allow_html=True)
            button3 = st.button('Make docking', on_click = callback3())
            if (button3 or st.session_state.button_clicked3):
                st.subheader('Protein-Ligand docking results using DiffDock solution')
                # with open("7wt1_complex.pdb", "r") as f:
                #     pdb_str = f.read()
                #
                # render_mol((pdb_str))
                image_orig = Image.open('images/orig_docked.png')
                image_mut8 = Image.open('images/mut8_docked.png')
                st.image(image_orig, caption='Docking for ORIGINAL protein')
                st.image(image_mut8, caption='Docking for MUTATED protein')



if __name__ == "__main__":
    main()
