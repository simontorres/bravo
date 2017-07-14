# todos los archvivos de python deben llevar esto al principio
# esto es para que el codigo sea compatible con futuras versiones de python
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# la maxima longitud de las lineas es 80 caracteres

# los comentarios de una linea deben ser como este. Empieza con gato (#) y un
# espacio. crtl + / comenta una linea

"""
los comentarios multilinea pueden ser como el comentario anterior, osea
cada linea empieza con gato y un espacio en blanco o como este
"""

# las funciones deben tener la siguiente forma


def mi_funcion(arg1, arg2, arg3=False):
    """Esta es una docstring

    Aca va una descripcion mas larga de lo que hace la funcion. Hay muchos mas
    atributos que se pueden incluir pero los siguientes son de los mas
    importantes

    Args:
        arg1 (str): Como describir un argumento tipo string
        arg2 (float): para un float
        arg3 (bool): Un argumento booleano

    Returns:
        Una descripcion del retorno

    """
    print('Hola Mundo')
    return True


# clases usan CamelCase, cada primera letra de cada palabra es mayuscula.
# ademas se recomienda usar palabras claras no cripticas y no muy larga tampoco


class MiClase(object):

    def __init__(self):
        self.name = 'Bravo'

    def saludar(self):
        print('Grande {:s}!'.format(self.name))