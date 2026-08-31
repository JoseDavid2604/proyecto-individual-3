# Proyecto Corto #3 - Control Automático

## Tecnológico de Costa Rica
### Escuela de Ingeniería Electrónica
### EL-5409 Laboratorio de Control Automático

**Estudiante:** Jose David Luna Herrera  
**Carné:** 2020114728  
**Profesor:** Ing. Luis C. Rosales  
**Semestre:** II Semestre 2026

---

## Descripción

En este proyecto corto se desarrolla un programa en MATLAB para el análisis y diseño de sistemas de control utilizando el método del **lugar de las raíces (Root Locus)**.

El programa permite ingresar los polos y ceros de una función de transferencia y, a partir de estos datos, construir la función de transferencia de la planta y obtener su ecuación característica.

Posteriormente, se muestra el lugar de las raíces y se permite seleccionar mediante click la ubicación deseada de los polos del sistema.

A partir de los polos seleccionados, el programa calcula un compensador dinámico y verifica que los polos obtenidos del sistema compensado coincidan con los polos deseados.

---

## Contenido

- Ingresa los polos y ceros de una función de transferencia.
- Construye la función de transferencia de la planta.
- Obtiene la ecuación característica del sistema.
- Grafica el lugar de las raíces.
- Permite seleccionar la ubicación deseada de los polos.
- Diseña un compensador dinámico a partir de los polos deseados.
- Obtiene la nueva ecuación característica del sistema compensado.
- Verifica la ubicación de los polos obtenidos.
- Compueba la respuesta al escalón de la planta y del sistema compensado.

---

## Funcionamiento del programa

Al ejecutar el programa en mathlab, se solicitan los siguientes datos:

1. **Ceros de la planta**
2. **Polos de la planta**
3. **Ganancia de la planta**

Por ejemplo:

```matlab
Ceros = []
Polos = [-1 -3 -5]
K = 1
