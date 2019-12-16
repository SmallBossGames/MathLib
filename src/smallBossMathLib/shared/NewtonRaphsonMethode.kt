package smallBossMathLib.shared

import kotlin.math.abs
import kotlin.math.sign

fun newtonRaphson(
    guess: Double,
    min: Double,
    max: Double,
    digits: Int,
    maxIter: Int,
    f:(value: Double) -> Pair<Double, Double>
) : Double {
    require(max > min)

    var min = min
    var max = max
    var guess = guess

    var f0 = 0.0;
    var f1 = 0.0;
    var result = guess

    val factor = (1.0 / (1 shl (digits - 1)))
    var delta = Double.MAX_VALUE
    var delta1 = Double.MAX_VALUE

    var maxRangeF = 0.0
    var minRangeF = 0.0
    var count = maxIter

    do {
        val delta2 = delta1
        delta1 = delta

        val functionResult = f(result)
        f0 = functionResult.first
        f1 = functionResult.second

        --count;

        if(f0 == 0.0)
            break

        if(f1 == 0.0)
            throw Exception("Derivative cannot be equal zero")
        else
            delta = f0/f1


        if(abs(delta * 2.0) > abs(delta2))
        {
            val shift = if (delta > 0.0)
                (result - min) / 2.0
            else
                (result - max) / 2.0

            delta = if((result != 0.0) && (abs(shift) > abs(result))) {
                sign(delta) * abs(result) * 1.1
            } else {
                shift
            }

            delta1 = 3 * delta
        }

        guess = result
        result -= delta

        if(result <= min)
        {
            delta = 0.5 * (guess - min)
            result = guess - delta
            if(result == min || result == max)
                break
        }
        else if (result >= max)
        {
            delta = 0.5 * (guess - max)
            result = guess - delta
            if(result == min || result == max)
                break
        }

        if (delta > 0)
        {
            max = guess
            maxRangeF = f0
        }
        else
        {
            min = guess
            minRangeF = f0
        }

        if(maxRangeF * minRangeF > 0)
            throw Exception("There appears to be no root to be found")

    } while (count != 0 && (abs(result * factor) < abs(delta)))

    return result
}