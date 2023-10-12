from Bio import Align
from Bio.Align import substitution_matrices
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "tu_correo_electronico@ejemplo.com"

# Número de acceso
accession_number1 = input("\nNumero de Secuencia 1: ")
accession_number2 = input("\nNumero de Secuencia 2: ")
print("\n")

# búsqueda en GenBank
handle = Entrez.efetch(
    db="nucleotide", id=accession_number1, rettype="gb", retmode="text"
)
# Lee la secuencia
record = SeqIO.read(handle, "genbank")
print("Número de acceso 1:", record.id)
print("Descripción 1:", record.description)
print("Longitud de la secuencia 1:", len(record.seq))
# print("Secuencia:")
# print(record.seq)

id1 = record.id
descrip1 = record.description
length1 = len(record.seq)
seq1 = record.seq

handle = Entrez.efetch(
    db="nucleotide", id=accession_number2, rettype="gb", retmode="text"
)

record = SeqIO.read(handle, "genbank")
print("Número de acceso 2:", record.id)
print("Descripción 2:", record.description)
print("Longitud de la secuencia 2:", len(record.seq))
# print("Secuencia:")
# print(record.seq)
id2 = record.id
descrip2 = record.description
length2 = len(record.seq)
seq2 = record.seq


# matriz BLOSUM62 y PAM250 se cambia desde el bloque
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

# Realiza el alineamiento global
global_alignments = aligner.align(seq1, seq2)

# Configura el alineador para realizar un alineamiento local
aligner.mode = "local"

# Realiza el alineamiento local
local_alignments = aligner.align(seq1, seq2)

i = 0
# Imprime los alineamientos globales
print("Alineamientos Globales:")
for alignment in global_alignments:
    print(alignment)

    global_identity = (
        sum(a == b for a, b in zip(alignment[0], alignment[1])) / len(alignment[0])
    ) * 100
    global_similarity = (
        alignment.score / max(len(alignment[0]), len(alignment[1]))
    ) * 100
    print(f"Score Total Cadena (Global): {alignment.score:.2f}")
    print(f"Porcentaje de identidad (Global): {global_identity:.2f}%")
    print(f"Porcentaje de similitud (Global): {global_similarity:.2f}%")
    if i < max(len(alignment[0]), len(alignment[1])):
        break
# Imprime los alineamientos locales
print("\nAlineamientos Locales:")
for alignment in local_alignments:
    print(alignment)

    local_identity = (
        sum(a == b for a, b in zip(alignment[0], alignment[1])) / len(alignment[0])
    ) * 100
    local_similarity = (
        alignment.score / max(len(alignment[0]), len(alignment[1]))
    ) * 100
    print(f"Score Total Cadena (Global): {alignment.score:.2f}")
    print(f"Porcentaje de identidad (Local): {local_identity:.2f}%")
    print(f"Porcentaje de similitud (Local): {local_similarity:.2f}%")
    if i < max(len(alignment[0]), len(alignment[1])):
        break
