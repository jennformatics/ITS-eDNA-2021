from matplotlib_venn import venn3
import matplotlib.pyplot as plt

f = open('presences.txt')
l = [x.strip().split() for x in f.readlines()]
f.close()

fulldict = {}
for m in l:
  s = m[6]
  fulldict[s] = {}
  
  m1 = m[1].strip(';').split(';')
  if len(m1) > 1:
    fulldict[s][m[0]] = (m1[0], int(m1[1][5:]))
  else:
    fulldict[s][m[0]] = ('', 0)
  
  m3 = m[3].strip(';').split(';')
  if len(m3) > 1:
    fulldict[s][m[2]] = (m3[0], int(m3[1][5:]))
  else:
    fulldict[s][m[2]] = ('', 0)
  
  m5 = m[5].strip(';').split(';')
  if len(m5) > 1:
    fulldict[s][m[4]] = (m5[0], int(m5[1][5:]))
  else:
    fulldict[s][m[4]] = ('', 0)

otusets = fulldict[fulldict.keys()[10]].keys()

def regenAndVenn(min=0):
  setdict = dict(zip(otusets, [[], [], []]))
  
  for s in otusets:
    for seq in fulldict:
      if fulldict[seq][s][1] > min:
        setdict[s].append(seq)
    setdict[s] = set(setdict[s])
  
  plt.figure(figsize=(4,4))
  v = venn3(setdict.values(), set_labels = setdict.keys())
  plt.title("Sample Venn diagram")
  plt.show()



def get_venn_sections(sets):
    """
    From https://stackoverflow.com/a/42161940/65664

    Given a list of sets, return a new list of sets with all the possible
    mutually exclusive overlapping combinations of those sets.  Another way
    to think of this is the mutually exclusive sections of a venn diagram
    of the sets.  If the original list has N sets, the returned list will
    have (2**N)-1 sets.
    
    Parameters
    ----------
    sets : list of set
    
    Returns
    -------
    combinations : list of tuple
        tag : str
            Binary string representing which sets are included / excluded in
            the combination.
        set : set
            The set formed by the overlapping input sets.
    """
    num_combinations = 2 ** len(sets)
    bit_flags = [2 ** n for n in range(len(sets))]
    flags_zip_sets = [z for z in zip(bit_flags, sets)]
    
    combo_sets = []
    for bits in range(num_combinations - 1, 0, -1):
        include_sets = [s for flag, s in flags_zip_sets if bits & flag]
        exclude_sets = [s for flag, s in flags_zip_sets if not bits & flag]
        combo = set.intersection(*include_sets)
        combo = set.difference(combo, *exclude_sets)
        tag = ''.join([str(int((bits & flag) > 0)) for flag in bit_flags])
        combo_sets.append((tag, combo))
    return combo_sets