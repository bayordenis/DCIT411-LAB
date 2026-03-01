import matplotlib.pyplot as plt

names = ["Needleman", "Smith-Waterman"]
values = [78.4, 54.6]  # replace with real numbers from pairwise.py

plt.bar(names, values)
plt.ylabel("Identity %")
plt.title("Pairwise Identity Comparison")
plt.savefig("results/figures/identity_comparison.png")