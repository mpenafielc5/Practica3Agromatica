# SESIÓN 3: Visualización con Matplotlib
# --------------------------------------

from Bio import SeqIO
import matplotlib.pyplot as plt

ids = []
longitudes = []

for record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    ids.append(record.id)
    longitudes.append(len(record.seq))

# Gráfico de barras (primeras 10 secuencias)
plt.bar(ids[:10], longitudes[:10])
plt.xticks(rotation=90)
plt.ylabel("Longitud")
plt.title("Longitud de secuencias (primeras 10)")
plt.tight_layout()
plt.savefig("grafico_longitudes.png")
plt.show()

# Gráfico de pastel de composición A, T, G, C (primera secuencia)
seq = str(list(SeqIO.parse("ls_orchid.fasta.txt", "fasta"))[0].seq)
conteo = [seq.count("A"), seq.count("T"), seq.count("G"), seq.count("C")]
bases = ["A", "T", "G", "C"]

plt.pie(conteo, labels=bases, autopct="%1.1f%%", startangle=90)
plt.title("Composición de bases (primera secuencia)")
plt.savefig("grafico_bases.png")
plt.show()
