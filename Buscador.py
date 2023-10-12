from Bio import Entrez
from Bio import SeqIO

Entrez.email = "tu_correo_electronico@ejemplo.com"

# Número de acceso
accession_number = input("\nNumero de Secuencia 1: ")

# búsqueda en GenBank
handle = Entrez.efetch(
    db="nucleotide", id=accession_number, rettype="gb", retmode="text"
)

# Lee la secuencia
record = SeqIO.read(handle, "genbank")

# Imprime información sobre la secuencia
print("Número de acceso:", record.id)
print("Descripción:", record.description)
print("Longitud de la secuencia:", len(record.seq))
print("Secuencia:")
print(record.seq)
