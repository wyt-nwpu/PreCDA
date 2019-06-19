package javasrc;

import java.io.Serializable;
import java.math.BigDecimal;

public class DoubleUtil implements Serializable {
	private static final long serialVersionUID = -3345205828566485102L;
	private static final Integer DEF_DIV_SCALE = 2;

	public static Double add(Number value1, Number value2) {
		BigDecimal b1 = new BigDecimal(Double.toString(value1.doubleValue()));
		BigDecimal b2 = new BigDecimal(Double.toString(value2.doubleValue()));
		return b1.add(b2).doubleValue();
	}

	public static double sub(Number value1, Number value2) {
		BigDecimal b1 = new BigDecimal(Double.toString(value1.doubleValue()));
		BigDecimal b2 = new BigDecimal(Double.toString(value2.doubleValue()));
		return b1.subtract(b2).doubleValue();
	}

	public static Double mul(Number value1, Number value2) {
		BigDecimal b1 = new BigDecimal(Double.toString(value1.doubleValue()));
		BigDecimal b2 = new BigDecimal(Double.toString(value2.doubleValue()));
		return b1.multiply(b2).doubleValue();
	}
	
	public static Double mul(String value1, String value2) {
		BigDecimal b1 = new BigDecimal(value1);
		BigDecimal b2 = new BigDecimal(value2);
		return b1.multiply(b2).doubleValue();
	}

	public static Double div(Double dividend, Double divisor) {
		return div(dividend, divisor, DEF_DIV_SCALE);
	}

	public static Double div(Double dividend, Double divisor, Integer scale) {
		if (scale < 0) {
			throw new IllegalArgumentException(
					"The scale must be a positive integer or zero");
		}
		BigDecimal b1 = new BigDecimal(Double.toString(dividend));
		BigDecimal b2 = new BigDecimal(Double.toString(divisor));
		return b1.divide(b2, scale, BigDecimal.ROUND_HALF_UP).doubleValue();
	}

	public static Double round(Double value, Integer scale) {
		if (scale < 0) {
			throw new IllegalArgumentException(
					"The scale must be a positive integer or zero");
		}
		BigDecimal b = new BigDecimal(Double.toString(value));
		BigDecimal one = new BigDecimal("1");
		return b.divide(one, scale, BigDecimal.ROUND_HALF_UP).doubleValue();
	}
}
