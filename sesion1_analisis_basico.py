# SESIÓN 1: Análisis básico de secuencias
# ---------------------------------------

# 1. Abrir y leer el archivo FASTA
archivo = open("ls_orchid.fasta.txt", "r")
contenido = archivo.read()
print(contenido[:500])  # Solo muestra los primeros 500 caracteres
archivo.close()

# 2. Copiar una secuencia manualmente (tomada del archivo)
secuencia = "ATGCGTACGTAGCTAGCTAGCTA"

# 3. Funciones simples de análisis

def longitud(seq):
    return len(seq)

def conteo_bases(seq):
    return {
        "A": seq.count("A"),
        "T": seq.count("T"),
        "C": seq.count("C"),
        "G": seq.count("G")
    }

def porcentaje_gc(seq):
    g = seq.count("G")
    c = seq.count("C")
    return round(((g + c) / len(seq)) * 100, 2)

# 4. Mostrar resultados
print("\n--- RESULTADOS DEL ANALISIS ---")
print("Longitud de la secuencia:", longitud(secuencia))
print("Conteo de bases:", conteo_bases(secuencia))
print("Porcentaje GC:", porcentaje_gc(secuencia), "%")
