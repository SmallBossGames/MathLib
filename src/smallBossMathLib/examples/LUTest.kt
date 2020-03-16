package smallBossMathLib.examples

import smallBossMathLib.shared.Matrix2D

fun luTest(){
    val matrix = Matrix2D(4)
    matrix[0,0] = 1.0
    matrix[0, 1] = 2.0
    matrix[0, 2] = 4.0
    matrix[0, 3] = 1.0
    matrix[1,0] = 2.0
    matrix[1, 1] = 8.0
    matrix[1, 2] = 6.0
    matrix[1, 3] = 4.0
    matrix[2,0] = 3.0
    matrix[2, 1] = 10.0
    matrix[2, 2] = 8.0
    matrix[2, 3] = 8.0
    matrix[3,0] = 4.0
    matrix[3, 1] = 12.0
    matrix[3, 2] = 10.0
    matrix[3, 3] = 6.0

    val rightPart = doubleArrayOf(21.0, 52.0, 79.0, 82.0)
    val result = DoubleArray(rightPart.size)

    matrix.makeLU()
    matrix.solveLU(rightPart, result)

    val inverseMatrix = Matrix2D(4)
    matrix.inverseLU(inverseMatrix)

    println(matrix.toString())
    println("Det: ${matrix.detLU()}")
    println("Inverse: $inverseMatrix")

    for (i in result)
        println(i)
}