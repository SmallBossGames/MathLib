package smallBossMathLib.shared

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

    val innerArray = Array(size) {DoubleArray(size)}

    inline val indices get() = innerArray.indices

    operator fun get(x: Int, y: Int) = innerArray[y][x]
    operator fun set(x: Int, y: Int, value: Double) {
        innerArray[y][x] = value
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
        for (p in 0 until size) {
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

        for (i in rightPart.indices){
            outResultVector[i] = rightPart[i]
        }

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

        for (col in outInverseMatrix.innerArray){
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