import matplotlib.pyplot as plt
from polytrope import polytrope

for n in [1.5, 3]:
    p = polytrope(n)
    plt.plot(p.z, p.w, label='n='+str(n))

plt.axis([0,7,0,1.05])
plt.xlabel('z', fontsize=14)
plt.ylabel('w(z) = solution of equation (9)', fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
