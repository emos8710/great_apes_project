import matplotlib.pyplot as plt

pos = range(1, 51)
cov = (60, 80, 95, 100, 95, 80, 70, 50, 30, 25, 15, 15, 15, 15, 15, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 5, 5, 5,
	   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 20, 30, 40, 35, 30, 20)

print pos
print len(pos), len(cov)

plt.figure()
plt.bar(pos, cov)
plt.show()
