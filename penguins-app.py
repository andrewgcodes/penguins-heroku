import streamlit as st
import streamlit.components.v1 as components
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.use("Agg")
from Bio.Seq import Seq 
from Bio import SeqIO
from collections import Counter
import neatbio.sequtils as utils
import numpy as np 
from PIL import Image 
import requests as req
from stmol import component_3dmol
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import kinetics
import pubchempy as pcp
from pubchempy import get_compounds, Compound

def lipinski(smile):
	# Convert into Chem object
	mol = Chem.MolFromSmiles(smile)

	MolWt = Descriptors.MolWt(mol)
	MolLogP = Descriptors.MolLogP(mol)
	NumHDonors = Lipinski.NumHDonors(mol)
	NumHAcceptors = Lipinski.NumHAcceptors(mol)

	return NumHDonors, NumHAcceptors, MolWt, MolLogP

def delta(x,y):
    return 0 if x == y else 1


def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]


def plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)


def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice


# Convert to Fxn
def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # on y-axis list all sequences of seq 1
    yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()


def gc_content(seq):
    result = float(str(seq).count('G') + str(seq).count('C'))/len(seq) * 100
    return result

def at_content(seq):
    result = float(str(seq).count('A') + str(seq).count('T'))/len(seq) * 100
    return result



def main():
    """A Simple Streamlit App """
    st.title("BioInformatics App")
    st.set_option('deprecation.showfileUploaderEncoding', False)

    activity = ['Intro','SequenceAnalysis','DotPlot','ProteinSearch',"MoleculeVisualizer", "ChemicalSearch"]
    choice = st.sidebar.selectbox("Select Activity",activity)
    if choice == 'Intro':
        st.subheader("Intro")
        st.write(""" This is a bioinformatics web app made with Python and Streamlit. Use the left panel dropdown to choose the various features to use. """)
    elif choice == "SequenceAnalysis":
        st.subheader("DNA Sequence Analysis")

        seq_file = st.file_uploader("Upload FASTA File",type=["fasta","fa"])

        if seq_file is not None:
            dna_record = SeqIO.read(seq_file,"fasta")
            # st.write(dna_record)
            dna_seq = dna_record.seq

            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record.description)
            elif details == "Sequence":
                st.write(dna_record.seq)


            # Nucleotide Frequencies
            st.subheader("Nucleotide Frequency")
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)
            adenine_color = st.beta_color_picker("Adenine Color")
            thymine_color = st.beta_color_picker("thymine Color")
            guanine_color = st.beta_color_picker("Guanine Color")
            cytosil_color = st.beta_color_picker("cytosil Color")

            if st.button("Plot Freq"):
                barlist = plt.bar(dna_freq.keys(),dna_freq.values())
                barlist[2].set_color(adenine_color)
                barlist[3].set_color(thymine_color)
                barlist[1].set_color(guanine_color)
                barlist[0].set_color(cytosil_color)

                st.pyplot()

            st.subheader("DNA Composition")
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.json({"GC Content":gc_score,"AT Content":at_score})

            # Nucleotide Count
            nt_count = st.text_input("Enter Nucleotide Here","Type Nucleotide Alphabet")
            st.write("Number of {} Nucleotide is ::{}".format((nt_count),str(dna_seq).count(nt_count)))

            # Protein Synthesis
            st.subheader("Protein Synthesis")
            p1 = dna_seq.translate()
            aa_freq = Counter(str(p1))

            if st.checkbox("Transcription"):
                st.write(dna_seq.transcribe())

            elif st.checkbox("Translation"):
                st.write(dna_seq.translate())

            elif st.checkbox("Complement"):
                st.write(dna_seq.complement())

            elif st.checkbox("AA Frequency"):
                st.write(aa_freq)

            elif st.checkbox("Plot AA Frequency"):
                aa_color = st.beta_color_picker("Pick An Amino Acid Color")
                # barlist = plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                # barlist[2].set_color(aa_color)
                plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                st.pyplot()

            elif st.checkbox("Full Amino Acid Name"):
                aa_name = str(p1).replace("*","")
                aa3 = utils.convert_1to3(aa_name)
                st.write(aa_name)
                st.write("=====================")
                st.write(aa3)

                st.write("=====================")
                st.write(utils.get_acid_name(aa3))
                
    elif choice == "ProteinSearch":
        st.subheader("Search for Papers Related to a Protein")
        st.write(""" Try entering ACE2 and coronavirus!""")

        ace2 = st.text_input("Query Protein")
        disease = st.text_input("Query Specifier (more specific thing to narrow down papers with)")

        if ace2 and disease is not None:
            protein = req.get('https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=10&gene='+ace2+'&organism=homo%20sapiens', headers = {'Accept':"application/json"})
            for i,v in enumerate(protein.json()[0]['references']):
                counter = 1
                try:
                    title = protein.json()[0]['references'][i]['citation']['title']
                    if counter ==10:
                        break
            
                    if title.find(disease) != -1:
                        st.write(title)
                        counter +=1
                except:
                    pass
            
    elif choice == "DotPlot":
        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file1 = st.file_uploader("Upload 1st FASTA File",type=["fasta","fa"])
        seq_file2 = st.file_uploader("Upload 2nd FASTA File",type=["fasta","fa"])

        if seq_file1 and seq_file2 is not None:
            dna_record1 = SeqIO.read(seq_file1,"fasta")
            dna_record2 = SeqIO.read(seq_file2,"fasta")
            # st.write(dna_record)
            dna_seq1 = dna_record1.seq
            dna_seq2 = dna_record2.seq

            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record1.description)
                st.write("=====================")
                st.write(dna_record2.description)
            elif details == "Sequence":
                st.write(dna_record1.seq)
                st.write("=====================")
                st.write(dna_record2.seq)


            cus_limit = st.number_input("Select Max number of Nucleotide",10,200,50)
            if st.button("Dot Plot"):
                st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                dotplotx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit])

                st.pyplot()



        



    elif choice == "MoleculeVisualizer":
        st.subheader("Look at a molecule! Pre-loaded example is the Covid-19 Spike Protein.")

        component_3dmol()
    
    elif choice == "ChemicalSearch":
        st.title("Search for chemicals and get info. Pre-loaded example: imatinib")
        user_compound = st.text_input("Enter compound name", 'imatinib')
        if user_compound is not None:
            results = pcp.get_compounds(user_compound,'name')
            for compound in results:
                st.write('Compound ID: '+ str(compound.cid))
                st.write('SMILES: '+compound.isomeric_smiles)
        
                vioxx = Compound.from_cid(compound.cid)
                st.write('Molecular Formula: '+vioxx.molecular_formula)
                st.write('Molecular Weight: '+str(vioxx.molecular_weight))
                st.write('IUPAC Name: '+vioxx.iupac_name)
                st.write('xlogp value: '+str(vioxx.xlogp))
                
        st.header("Input SMILE")
        user_smile = st.text_input("Enter text below")

        # Calculation
        hDonar = (lipinski(user_smile)[0])
        hAccep = (lipinski(user_smile)[1])
        molWgt = (lipinski(user_smile)[2])
        logPVa = (lipinski(user_smile)[3])

        st.header("Lipinski's Descriptors Values")

        st.write(pd.DataFrame({

            'H Donars': pd.Series(hDonar),
            'H Acceptors': pd.Series(hAccep),
            'Molecular Mass (Dalton)': pd.Series(molWgt),
            'LogP': pd.Series(logPVa)
                }))


if __name__ == '__main__':
    main()
