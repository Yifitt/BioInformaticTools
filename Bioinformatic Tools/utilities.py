from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd
import time

def read_csv_to_dataframe(filepath):
    try:
        df = pd.read_csv(filepath)
        print(f"Successfully read '{filepath}' into a DataFrame.")
        return df
    except FileNotFoundError:
        print(f"Error: The file at '{filepath}' was not found. Please check the path.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading the CSV: {e}")
        return None

def readTextFile(filePath):
    with open(filePath, 'r') as f:
        return "".join([l.strip() for l in f.readlines()])


def writeTextFile(filePath, seq, mode='w'):
    with open(filePath, mode) as f:
        f.write(seq + '\n')


def read_FASTA(filePath):
    with open(filePath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]

    FASTADict = {}
    FASTALabel = ""

    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line

    return FASTADict

def write_FASTA(filePath,seq_title ,seq, mode='w'):
    with open(filePath, mode) as f:
        f.write('>' + seq_title + '\n')
        f.write(seq + '\n')

def sendToBlast(filePath,delay=1):
    query = SeqIO.read(filePath,format="fasta")
    blast_file_name = "my_blast.xml"
    success = False
    while not success:
        try:
            result_handle = NCBIWWW.qblast("blastn", "nt", query.seq)

            blast_file = open(blast_file_name, "w")
            blast_file.write(result_handle.read())
            result_handle = open(blast_file_name)

            if (result_handle):
                print("XML file created")
                success = True
                return blast_file_name
        except Exception as e:
            print(f"Error: {e}")
            time.sleep(delay)

        finally:
            blast_file.close()
            result_handle.close()

    return None


class BlastProperties():
    def __init__(self,xml_path):
        self.results = []
        with open(xml_path) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        record = {
                            "title": alignment.title,
                            "score": hsp.score,
                            "e_value": hsp.expect,
                            "sequence": hsp.sbjct,
                            "matches": hsp.identities,
                            "alignment_length": hsp.align_length,
                            "gaps": hsp.gaps,
                            "taxonomy": alignment.title.split()[1] if len(alignment.title.split()) > 1 else "Unknown"
                        }
                        self.results.append(record)

    @property
    def all_results(self):
        if not self.results:
            print("Warning: No BLAST results to display in DataFrame.")
            return pd.DataFrame()
        return pd.DataFrame(self.results)




