#pragma once

class Color
{
public:
	float r, g, b;
	Color() : r(0), g(0), b(0) {}
	Color(float _r, float _g, float _b) : r(_r), g(_g), b(_b) {}
	Color(float _c) : r(_c), g(_c), b(_c) {}
	Color& operator + (const Color& c) { return Color(r + c.r, g + c.g, b + c.b); }
	Color& operator - (const Color& c) { return Color(r - c.r, g - c.g, b - c.b); }
	Color& operator * (const float d) { return Color(r * d, g * d, b * d); }
	const Color& operator * (const Color& c) { return Color(r * c.r, g * c.g, b * c.b); }
	Color& operator / (const Color& c) { return Color(r / c.r, g / c.g, b / c.b); }
	Color& operator / (const float d) { return Color(r / d, g / d, b / d); }
	const Color& operator += (const Color& c) { r += c.r; g += c.g; b += c.b; return *this; }
	bool operator > (const Color& c) { return (r > c.r && g > c.g && b > c.b); }
	bool operator < (const Color& c) { return (r < c.r && g < c.g && b < c.b); }
	void clamp(float minVal = 0.0, float maxVal = 1.0)
	{
		if (r < minVal)
			r = minVal;
		else if (r > maxVal)
			r = maxVal;

		if (g < minVal)
			g = minVal;
		else if (g > maxVal)
			g = maxVal;

		if (b < minVal)
			b = minVal;
		else if (b > maxVal)
			b = maxVal;
	}
};

template <class T>
T clamp(T val, T min = 0.0, T max = 1.0)
{
	if (val < min)
		val = min;
	else if (val > max)
		val = max;
	return val;
}

Color min(const Color a, const Color b)
{
	Color c;
	c.r = std::min(a.r, b.r);
	c.g = std::min(a.g, b.g);
	c.b = std::min(a.b, b.b);
	return c;
}

Color max(const Color a, const Color b)
{
	Color c;
	c.r = std::max(a.r, b.r);
	c.g = std::max(a.g, b.g);
	c.b = std::max(a.b, b.b);
	return c;
}

void RGBtoHSV(int r, int g, int b, double &h, double &s, double &v)
{

#define maximum(x, y, z) ((x) > (y)? ((x) > (z)? (x) : (z)) : ((y) > (z)? (y) : (z))) 
#define minimum(x, y, z) ((x) < (y)? ((x) < (z)? (x) : (z)) : ((y) < (z)? (y) : (z)))

	double red, green, blue;
	double max, min, delta;

	red = r / 255.0; green = g / 255.0; blue = b / 255.0;  /* r, g, b to 0 - 1 scale */

	max = maximum(red, green, blue);
	min = minimum(red, green, blue);

	v = max;        /* value is maximum of r, g, b */

	if (max == 0) {    /* saturation and hue 0 if value is 0 */
		s = 0;
		h = 0;
	}
	else {
		s = (max - min) / max;           /* saturation is color purity on scale 0 - 1 */

		delta = max - min;
		if (delta == 0)                    /* hue doesn't matter if saturation is 0 */
			h = 0;
		else {
			if (red == max)                  /* otherwise, determine hue on scale 0 - 360 */
				h = (green - blue) / delta;
			else if (green == max)
				h = 2.0 + (blue - red) / delta;
			else /* (blue == max) */
				h = 4.0 + (red - green) / delta;
			h = h * 60.0;
			if (h < 0)
				h = h + 360.0;
		}
	}
}

void HSVtoRGB(int& R, int& G, int& B, double h, double s, double v)
{
	double hh, p, q, t, ff;
	long i;
	double r, g, b;
	if (s <= 0.0) // < is bogus, just shuts up warnings
	{
		r = v;
		g = v;
		b = v;
		R = r * 255;
		G = g * 255;
		B = b * 255;
		return;
	}

	hh = h;
	if (hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = v * (1.0 - s);
	q = v * (1.0 - (s * ff));
	t = v * (1.0 - (s * (1.0 - ff)));

	switch (i)
	{
	case 0:
		r = v;
		g = t;
		b = p;
		break;
	case 1:
		r = q;
		g = v;
		b = p;
		break;
	case 2:
		r = p;
		g = v;
		b = t;
		break;
	case 3:
		r = p;
		g = q;
		b = v;
		break;
	case 4:
		r = t;
		g = p;
		b = v;
		break;
	case 5:
	default:
		r = v;
		g = p;
		b = q;
		break;
	}
	R = r * 255;
	G = g * 255;
	B = b * 255;
}

double randomGen(double min = 0.0, double max = 1.0)
{
	static std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}

double convex(double x, double y)
{
	double f = 2 * (x - 100) - y;
	double g = (x + 200) - 2 * y;
	double h = 400 - y;
	double k = 400 - x;
	return std::max(std::max(-f, g), std::max(-h, -k));
}

double function(double x, double y)
{
	return -(100 * sin(x / 50) + x - y);
}

double blobby(double x, double y)
{
	double f = 100 * sin(x / 50 + 400) + x - y;
	double g = 200 * cos(x / 100 - 70) + 250 - y;
	double h = 100 * sin(x / 200) - y;
	double k = 400 * cos(x / 400) - y;
	return std::max(std::max(-f, g), std::max(h, -k));
}

double star(double x, double y)
{
	double f = x + 4000 - 10 * y;
	double g = -4 * x + 4500 - 10 * y;
	double h = -2 * x + 1400 - y;
	double k = 2 * x - 200 - y;
	double m = 5 * x + 1000 - 10 * y;
	double t1 = std::max(std::max(-f, g), m);
	double t2 = std::max(std::max(-h, -k), m);
	double t3 = std::max(std::max(g, -h), -k);
	return std::min(std::min(t1, t2), t3);
}

Color shapeColor(double functionVal)
{
	Color c;
	if (functionVal <= 0.0)
	{
		c.r += 226.0 / 256.0;
		c.g += 48.0 / 256.0;
		c.b += 81.0 / 256.0;
	}
	else
	{
		c.r += 65.0 / 256.0;
		c.g += 76.0 / 256.0;
		c.b += 83.0 / 256.0;
	}
	return c;
}

Color shaded(double x, double y)
{
	Color c;
	double f = pow(x - 400, 2) + pow(y - 400, 2) - pow(200, 2);
	double g = pow(x - 350, 2) + pow(y - 350, 2) - pow(100, 2);
	if (f <= 0)
	{
		c.r = (0 - f / pow(200, 2)) * 226.0 / 256.0;
		c.g = (0 - f / pow(200, 2)) * 48.0 / 256.0;
		c.b = 81.0 / 256.0;
		if (g <= 0)
		{
			c.r += pow((0 - g / pow(100, 2)), 2);
			c.g += pow((0 - g / pow(100, 2)), 2);
			c.b += pow((0 - g / pow(100, 2)), 2);
		}
	}
	else
	{
		c.r = 65.0 / 256.0;
		c.g = 76.0 / 256.0;
		c.b = 83.0 / 256.0;
	}
	c.r = clamp<double>(c.r);
	c.g = clamp<double>(c.g);
	c.b = clamp<double>(c.b);

	return c;
}
