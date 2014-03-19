import brewer2mpl

# "primary" type colors with red, blue, etc
set1 = brewer2mpl.get_map('Set1', 'Qualitative', 9)

# Map color names to colors in set1:
mpl_to_set1 = {'red': set1[0],
               'blue': set1[1],
               'green': set1[2],
               'purple': set1[3],
               'orange': set1[4],
               'yellow': set1[5],
               'brown': set1[6],
               'pink': set1[7],
               'grey': set1[8],
               'gray': set1[8]}

# More colorblind friendly colors but harder to describe like "teal", "peach",
# "gray-blue." Lighter than the dark2 colors. Good for larger item plotting
# like barplots
set2 = brewer2mpl.get_map('Set2', 'Qualitative', 8)

# Darker versions of the same colors in set2. Good for small-sized plotting
# items like scatterplots
dark2 = brewer2mpl.get_map('Dark2', 'Qualitative', 8)

