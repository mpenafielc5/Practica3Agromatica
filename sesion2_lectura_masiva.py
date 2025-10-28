# SESIÓN 2: Lectura automatizada con Biopython
# --------------------------------------------

from Bio import SeqIO

longitudes = []

for record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print("ID:", record.id)
    print("Longitud:", len(record.seq))
    longitudes.append(len(record.seq))

# Estadísticas generales
num_secuencias = len(longitudes)
promedio = sum(longitudes) / num_secuencias
max_longitud = max(longitudes)

print("\n--- ESTADISTICAS GENERALES ---")
print("Numero total de secuencias:", num_secuencias)
print("Longitud promedio:", round(promedio, 2))
print("Secuencia mas larga:", max_longitud)
