Each file corresponds to a matrix of the same name.

Each line in the file corresponds to a character in the matrix.

The first letter of each line determines whether the character is:

- I: Inverse Neomorphic - coded as '0: present, 1: absent' instead of '0: absent, 1: present'
- N: Neomorphic
- T: Transformational
- X: Uncoded (as it doesn't contain any inapplicable characters)

Subsequent characters of the line are not parsed; I tend to specify the character number 
[in `[]`s] after a space, followed by any comments (remembering not to introduce line breaks).

Template files can be generated using the perl script `character_type_template.pl`
