package smallBossMathLib.implicitDifferentialEquations

import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.linear.RealMatrix
import smallBossMathLib.shared.createUnitMatrix
import smallBossMathLib.shared.findJacobiMatrix

const val a = 0.43586652150846

class MK22Integrator {

    fun integrate(n: Int){
        //val J = Matrix
        val e = createUnitMatrix(n)
        //val j = findJacobiMatrix()

        //e.po
        //val d = e.subtract()
    }


}