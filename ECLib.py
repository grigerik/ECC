# -*- coding: utf-8 -*-


class ECurve(object):
    # класс эллиптической кривой над полем простого числа
    # y**2 = x**3 + a*x + b % p

    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p
        self.discriminant = 4*self.a**3 + 27*self.b**2
        self.isSmooth = True if self.discriminant != 0 else False
        if not self.isSmooth:
            raise Exception("Эллиптическая кривая не гладкая")

    def test_point(self, x, y):
        return (y**2) % self.p == (x**3 + self.a*x+self.b) % self.p

    def __str__(self):
        a = str(self.a)
        b = str(self.b)
        p = str(self.p)
        if a != "0":
            if b != "0":
                return "y^2 = x^3 + %s*x + %s mod %s" % (a, b, p)
            return "y^2 = x^3 + %s*x mod %s" % (a, p)
        if b != "0":
            return "y^2 = x^3 + %s mod %s" % (b, p)
        return "y^2 = x^3 mod %s" % p

    def __eq__(self, other):
        return (self.a, self.b, self.p) == (other.a, other.b, other.p)


class Ideal(object):
    def __init__(self, curve):
        self.curve = curve

    def __str__(self):
        return "Идеальная точка"

    def __neg__(self):
        return self

    def __add__(self, other):
        if self.curve != other.curve:
            raise TypeError("Идеальная точка другой группы точек эллиптической кривой")
        return other

    def __mul__(self, n):
        if not isinstance(n, (int, long)):
            raise Exception("Множитель точки не является числовым типом")
        return self

    def __rmul__(self, n):
        return self * n

    def __eq__(self, other):
        if isinstance(other, Ideal) and self.curve == other.curve:
            return True
        return False


class Point(Ideal):

    def __init__(self, curve, x, y):
        super(Point, self).__init__(curve)
        self.x = x
        self.y = y
        if not curve.test_point(x, y):
            raise Exception("Точка (%s;%s) не принадлежит группе точек эллиптической кривой" % (x, y))

    def __str__(self):
        return "Точка (%s;%s) принадлежит группе точек эллиптической кривой %s" % (self.x, self.y, self.curve)

    def __neg__(self):
        return Point(self.curve, self.x, -self.y)

    def __add__(self, other):
        if type(other) is Ideal:
            return self
        if type(other) is not Point:
            raise TypeError("Слагаемое %s не точка" % other)
        if self.curve != other.curve:
            raise TypeError("Слагаемые не входят в одну и ту же группу эллиптической кривой")
        p = self.curve.p
        x1, y1, x2, y2 = self.x, self.y, other.x, other.y
        if x1 == x2:
            if y1 != y2 or y1 == 0:
                return Ideal(self.curve)
            m = ((3 * x1 ** 2 + self.curve.a) * reverse_a(p, 2*y1)) % p
        else:
            m = ((y2 - y1) * reverse_a(p, (x2 - x1))) % p
        x3 = (m * m - x2 - x1) % p
        y3 = (m * (x1 - x3) - y1) % p
        return Point(self.curve, x3, y3)

    def __mul__(self, n):
        if not isinstance(n, (int, long)):
            raise Exception("Множитель точки не является числовым типом")
        elif n == 0:
            return Ideal(self.curve)
        elif n < 0:
            return (-self)*(-n)
        else:
            result = Ideal(self.curve)
            addend = self
            for bit in bits(n):
                if bit == 1:
                    result += addend
                addend = addend + addend
            return result

    def __rmul__(self, n):
        return self * n

    def __eq__(self, other):
        return (self.curve, self.x, self.y) == (other.curve, other.x, other.y)


def bits(n):
    while n:
        yield n & 1
        n >>= 1


def eea(a, b):
    if a < b:
        a, b = b, a
    if not b:
        d, x, y = a, 1, 0
        return d, x, y
    x_2, x_1, y_2, y_1 = 1, 0, 0, 1
    while b:
        q = a/b
        r = a - b*q
        x = x_2 - q*x_1
        y = y_2 - q*y_1
        a, b, x_2, x_1, y_2, y_1 = b, r, x_1, x, y_1, y
    d, x, y = a, x_2, y_2
    return d, x, y


def reverse_a(prime, a):
    a = a % prime
    if a == 0:
        raise ZeroDivisionError("Нельзя найти обратный множитель по модулю для нуля")
    # так как p>a, то я поменял местами x,y на y,x
    gcd, y, x = eea(prime, a)
    assert gcd == 1
    assert (prime*y + a*x) == gcd
    if gcd != 1:
        return 0
    x = x % prime
    return x
