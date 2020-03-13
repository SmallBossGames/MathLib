package smallBossMathLib.examples

import smallBossMathLib.explicitDifferentialEquations.StabilityControlSecondOrderIntegrator
import java.io.File
import kotlin.math.*

private fun uIn1(t: Double) = 0.5 * sin(2000 * PI * t)
private fun uIn2(t: Double) = 2.0 * sin(20000 * PI * t)
private fun q(u:Double) : Double{
    val gamma = 40.67286402e-9
    val delta = 17.7493332
    return gamma * (E.pow(delta*u) - 1.0)
}
private fun mainFunction(t: Double, i: DoubleArray, o: DoubleArray){
    val C = 1.6e-8
    val Cs = 2e-12
    val Cp = 10e-8
    val Lh = 4.45
    val Ls1 = 0.002
    val Ls2 = 5e-4
    val Ls3 = 5e-4
    val R = 25000.0
    val Rp = 50.0
    val Rg1 = 36.3
    val Rg2 = 17.3
    val Rg3 = 17.3
    val Ri = 50.0
    val Rc = 600.0
    val Ud1 = i[2]-i[4]-i[6]-uIn2(t)
    val Ud2 = -i[3]+i[5]-i[6]- uIn2(t)
    val Ud3 = i[3]+i[4]+i[6]+ uIn2(t)
    val Ud4 = -i[2]-i[5]+i[6]+ uIn2(t)


    o[0] = (1.0/C) * (i[7]-0.5*i[9]+0.5*i[10]+i[13]-(1.0/R)*i[0])
    o[1] = (1.0/C)*(i[8]-0.5*i[11]+0.5*i[12]+i[14]-(1.0/R)*i[1])
    o[2] = (1.0/Cs)*(i[9]-q(Ud1)+q(Ud4))
    o[3] = (1.0/Cs)*(-i[10]+q(Ud2)-q(Ud3))
    o[4] = (1.0/Cs)*(i[11]+q(Ud1)+q(Ud3))
    o[5] = (1.0/Cs)*(i[12]+q(Ud2)+q(Ud4))
    o[6] = (1.0/Cp)*(-(1.0/Rp)*i[6]+q(Ud1)+q(Ud2)-q(Ud3)-q(Ud4))
    o[7] = -(1.0/Lh)*i[0]
    o[8] = -(1.0/Lh)*i[1]
    o[9] = (1.0/Ls2)*(0.5*i[0]-i[2]-Rg2*i[9])
    o[10] = (1.0/Ls3)*(-0.5*i[0]-i[3]-Rg3*i[10])
    o[11] = (1.0/Ls2)*(0.5*i[1]-i[4]-Rg2*i[11])
    o[12] = (1.0/Ls3)*(-0.5*i[1]-i[5]-Rg3*i[12])
    o[13] = (1.0/Ls1)*(-i[0]+uIn1(t)-(Ri+Rg1)*i[13])
    o[14] = (1.0/Ls1)*(-i[1]-(Rc+Rg1)*i[14])
}

fun ringModualtorExample(){
    val integrator = StabilityControlSecondOrderIntegrator( 10000, 0.000001)
    val builder = StringBuilder()
    val output = DoubleArray(15) {0.0}

    integrator.addStepHandler(){
        t, y -> builder.append("${t};${y[13]} \n");
    }

    integrator.enableEvaluationCountCheck(20000)

    integrator.integrate(0.0, output, 0.001, output, ::mainFunction)

    val writingText = builder.replace(Regex("[.]"), ",")

    File("data_modulator_test.csv ").writeText(writingText)
}