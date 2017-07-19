import asciitable


def cargar_datos(tabla):
    data = asciitable.read(tabla)
    jda = data['col1']
    maga = data['col2']
    erra = data['col3']
    return jda, maga, erra