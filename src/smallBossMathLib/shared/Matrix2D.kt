package smallBossMathLib.shared

import kotlin.math.abs

class Matrix2D(val size: Int) {
    companion object{
        @JvmStatic
        fun createUnitMatrix2D(size: Int) : Matrix2D
        {
            val e = Matrix2D(size)

            for (i in e.indices)
                e[i, i] = 1.0

            return e
        }
    }

    val columns = Array(size) {DoubleArray(size)}
    private val transpositionVector = IntArray(size)

    inline val indices get() = columns.indices

    operator fun get(x: Int, y: Int) = columns[y][x]
    operator fun set(x: Int, y: Int, value: Double) {
        columns[y][x] = value
    }

    operator fun plusAssign(other: Matrix2D){
        for (i in this.indices)
            for (j in this.indices)
                 this[i, j] += other[i, j]
    }

    operator fun minusAssign(other: Matrix2D){
        for (i in this.indices)
            for (j in this.indices)
                this[i, j] -= other[i, j]
    }

    operator fun timesAssign (constant: Double){
        for (i in this.indices)
            for (j in this.indices)
                this[i, j] *= constant
    }

    operator fun divAssign (constant: Double){
        for (i in this.indices)
            for (j in this.indices)
                this[i, j] /= constant
    }

    fun makeLU(){
        for (i in transpositionVector.indices)
            transpositionVector[i] = i

        for (p in 0 until size) {
            var maxIndex = p
            for (i in p+1 until size)
                if(abs(this[i, p]) > abs(this[maxIndex, p]))
                    maxIndex = i

            if(p != maxIndex){
                val temp = transpositionVector[p]
                transpositionVector[p] = transpositionVector[maxIndex]
                transpositionVector[maxIndex] = temp
                for (i in 0 until size){
                    val temp1 = this[p, i]
                    this[p, i] = this[maxIndex, i]
                    this[maxIndex, i] = temp1
                }
            }

            for (r in (p+1) until size){
                val m = this[r, p] / this[p, p]
                this[r, p] = m
                for(c in (p+1) until size) {
                    this[r, c] = this[r, c] - m * this[p, c]
                }
            }
        }
    }

    fun multiply(vector: DoubleArray, outVector: DoubleArray){
        if(outVector.size != vector.size || vector.size != size)
            throw IllegalArgumentException()

        for (i in outVector.indices){
            outVector[i] = 0.0
            for (j in outVector.indices){
                outVector[i] += vector[j] * this[i, j]
            }
        }
    }

    @Throws(IllegalArgumentException::class)
    fun solveLU(rightPart: DoubleArray, outResultVector: DoubleArray){
        if(rightPart.size != size || outResultVector.size != size) {
            throw IllegalArgumentException()
        }

        val temp = DoubleArray(outResultVector.size)

        for (i in rightPart.indices){
            temp[i] = rightPart[transpositionVector[i]]
        }

        for (i in rightPart.indices)
            outResultVector[i] = temp[i]

        for (i in indices){
            val temp = outResultVector[i]
            for (j in (i + 1) until outResultVector.size){
                outResultVector[j] -= temp * this[j, i]
            }
        }

        for (i in (rightPart.size-1) downTo  0){
            val temp = outResultVector[i] / this[i, i]
            outResultVector[i] = temp

            for (j in (i-1) downTo  0){
                outResultVector[j] -= temp * this[j, i]
            }
        }
    }

    fun detLU() : Double{
        var result = 1.0
        for (i in indices){
            result *= this[i, i]
        }
        return result
    }

    fun inverseLU(outInverseMatrix: Matrix2D){
        if(outInverseMatrix.size != size)
            throw IllegalArgumentException()

        for (i in outInverseMatrix.indices){
            for (j in outInverseMatrix.indices){
                outInverseMatrix[i, j] = if (i == j) 1.0 else 0.0
            }
        }

        for (col in outInverseMatrix.columns){
            solveLU(col, col)
        }
    }

    override fun toString(): String {
        val sb = StringBuilder()

        sb.append('[')
        for (i in indices)
        {
            sb.append('[')
            for (j in indices){
                sb.append(this[i, j].toString()).append(',')
            }
            sb.append(']')
        }
        sb.append(']')

        return sb.toString()
    }
}