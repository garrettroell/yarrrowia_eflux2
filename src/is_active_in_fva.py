# define a function to determine if a reaction is active
def is_active_in_fva(fva_row):
  return abs(fva_row.maximum) > 0.1 or abs(fva_row.minimum) > 0.1