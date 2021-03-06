def num2sig(num):
    """Convert a p-avalue to a signifiacne annotation.

    Significance annotations are defined as follows:

    | P-value              | Signifiacne |
    | -------------------- | ----------- |
    |   0.05 < P           | ns          |
    |   0.01 < P <= 0.05   | *           |
    |  0.001 < P <= 0.01   | **          |
    | 0.0001 < P <= 0.001  | ***         |
    |          P <= 0.0001 | ****        |

    Parameters
    ----------
    num : float
        P-value to be converted.

    Returns
    -------
    str
        Signifiance annotation.
    """
    if 0.05 < num:
        sig = 'ns'
    elif 0.01 < num <= 0.05:
        sig = '*'
    elif 0.001 < num <= 0.01:
        sig = '**'
    elif 0.0001 < num <= 0.001:
        sig = '***'
    else:
        sig = '****'
    return sig
