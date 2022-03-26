from polytrope import polytrope

p = None
for n in [0, 1, 1.5, 2, 3, 4]:
    p = polytrope(n, init=p)
    print("${0}$ & ${1:.5}$ & ${2:.5}$\\\\"
          .format(n, p.zn, -p.wp[-1]/p.zn))
