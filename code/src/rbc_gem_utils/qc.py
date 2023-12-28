from cobra.core.metabolite import element_re

def standardardize_metabolite_formulas(metabolite_formulas):
    """Standardize metabolite chemical formulas by reordering elements according to the Hill system.

    Parameters
    ----------
    metabolite_formulas : dict
        Contains identifiers as keys and strings representing chemical formulas as values.

    Returns
    -------
    updated_formulas : dict
        Contains original keys and standardized chemical formulas as values.

    Raises
    ------
    ValueError
        Occurs when metabolite(s) have missing or non-alphanumeric formulas.
    
    """
    invalid_formula = {
        key: formula 
        for key, formula in metabolite_formulas.items() 
        if not formula           # validate formula is neither `None` nor an empty string
        or not formula.isalnum() # validate formula is alphanumeric, enforces integer coefficients
    }
    if invalid_formula:
        raise ValueError(f"Missing or invalid formulas for the following: {list(invalid_formula)}")

    updated_formulas = {}
    for key, formula in metabolite_formulas.items():
        # Parse formula
        parsed = element_re.findall(formula)
        # Identify whether carbons are in the formula
        carbons = [(e, n) for e, n in parsed if e == "C"]
        if carbons:
            # Presence of carbons, organize first carbons, then hydrogens, and alphabetically sort the remaining.
            hydrogens = [(e, n) for e, n in parsed if e == "H"]
            other = [(e, n) for e, n in parsed if e not in {"C", "H"}]
            parsed = carbons + hydrogens + sorted(other)
        else:
            # Alphabetically sort, regardless of the presence of hydrogens.
            parsed = sorted(parsed)
        # Add formula to return.
        updated_formulas[key] = "".join(e if n == 1 else f"{e}{n}" for e, n in parsed)
    
    # Return updated formulas as a dictionary with original keys
    return updated_formulas

