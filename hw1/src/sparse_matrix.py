import numpy as np

from bisect import bisect_left

class RedkaMatrika:
    """ Sparce matrix data type
        A sparce matrix is represented with two arrays V and I, s.t. the following holds:
            V[i][j] = a_(i, I[i][j])

        V[i] holds the nonzero values from i-th row of matrix A.
        I[i] holds the column indices of the nonzero values from i-th row.
    """

    def __init__(self, V, I):
        self.V = []
        self.I = []
        n = len(V) # Square matrices

        # Validate the V and I array. By saying strict=True zip will also make sure the length match
        for row_values, row_indices in zip(V, I, strict=True):
            entries = {}
            for value, column in zip(row_values, row_indices, strict=True):
                column = column
                if not 0 <= column < n:
                    raise IndexError("Invalid row matrix")
                if column in entries:
                    raise ValueError("Same nonzero element saved twice")
                if value != 0: # In case a rouge zero element is in V, we can fix this easily
                    entries[column] = float(value)

            # Sort by index
            columns = sorted(entries)
            self.I.append(columns)
            # Values are also sorted by their position (index)
            self.V.append([entries[column] for column in columns])


    @property
    def shape(self):
        """Returns the size of the matrix"""

        n = len(self.V)
        return n, n


    def _check_position(self, row, column):
        n, _ = self.shape
        if not (0 <= row < n and 0 <= column < n):
            raise IndexError("Index is out of bounds")


    def __getitem__(self, position):
        """Returns an element: value = A[i,j]"""

        row, column = position
        self._check_position(row, column)

        # Since we need to actually find the matching column index stored in I which is sorted
        # Thus use binary search
        offset = bisect_left(self.I[row], column)
        # Validate the index
        if offset < len(self.I[row]) and self.I[row][offset] == column:
            return self.V[row][offset]
        return 0.0 # We are returning a zero value from the matrix


    def __setitem__(self, position, value):
        """Stores an element: A[i,j] = value"""

        row, column = position
        self._check_position(row, column)

        value = float(value)
        # Same logic as in __getitem__
        offset = bisect_left(self.I[row], column)
        present = offset < len(self.I[row]) and self.I[row][offset] == column

        # A nonzero element is now set to 0, so remove it and its corresponding index from I
        if present and value == 0:
            del self.I[row][offset]
            del self.V[row][offset]
        elif present:
            # Set the value
            self.V[row][offset] = value
        elif value != 0:
            # A zero element is now nonzero thus added it to the matrix representation
            # Using bisect_left means that we can use the offset for inserting here
            self.I[row].insert(offset, column)
            self.V[row].insert(offset, value)

    def __matmul__(self, vector):
        """For multiplying matrix A with a vector on the right side. y = A @ x"""

        x = np.asarray(vector, dtype=float)
        n, _ = self.shape
        if x.ndim != 1:
            raise ValueError("The vector must be of dimension 1")
        # NxN dot N=> N vector
        result = np.zeros(n, dtype=float)

        for row, (values, columns) in enumerate(zip(self.V, self.I, strict=True)):
            # We don't even have to operate with the zero elements since they contribute nothing (i.e., zero)
            # This is (Ax)_i = Sum_j V_i[j] x_(I[i,j])
            result[row] = np.dot(values, x[columns])
        return result