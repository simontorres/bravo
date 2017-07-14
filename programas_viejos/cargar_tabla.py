import asciitable


def cargar_datos(tabla):
    data = asciitable.read(tabla)
    jdA=data['col1']
    magA=data['col2']
    errA=data['col3']