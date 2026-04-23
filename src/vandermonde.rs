use lin_algebra::gf2_matrix::GF2Matrix;
use lin_algebra::matrix::MatrixTrait;
use std::ops::{Deref, DerefMut};

#[derive(Clone)]
pub struct VanderMonde {
    pub matrix: GF2Matrix,
}

impl From<GF2Matrix> for VanderMonde {
    fn from(matrix: GF2Matrix) -> Self {
        Self { matrix }
    }
}

impl Deref for VanderMonde {
    type Target = GF2Matrix;
    fn deref(&self) -> &Self::Target {
        &self.matrix
    }
}

impl DerefMut for VanderMonde {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.matrix
    }
}

impl VanderMonde {
    /// Constructs a new `VanderMonde` matrix from a given set of binary row vectors.
    ///
    /// # Arguments
    ///
    /// * `elements` - A vector of vectors, where each inner vector represents a row over GF(2).
    ///   Each element must be `0` or `1`, representing a Boolean value.
    ///
    /// # Returns
    ///
    /// A `VanderMonde` instance containing the underlying GF(2) matrix.
    ///
    /// # Example
    ///
    /// ```
    /// use algebraic_immunity::vandermonde::VanderMonde;
    ///
    /// let rows = vec![
    ///     vec![1, 0, 1],
    ///     vec![0, 1, 0],
    /// ];
    /// let vm = VanderMonde::new(rows);
    /// assert_eq!(vm.elements, vec![
    ///     vec![1, 0, 1],
    ///     vec![0, 1, 0],
    /// ]);
    /// ```
    pub fn new(elements: Vec<Vec<u8>>) -> Self {
        VanderMonde {
            matrix: GF2Matrix::new(elements),
        }
    }

    pub fn __repr__(&self) -> String {
        let rows: Vec<String> = self
            .elements
            .iter()
            .map(|row| format!("{:?}", row))
            .collect();
        format!("[{}]", rows.join(", "))
    }

    /// Add a row to a VanderMonde matrix.
    ///
    /// # Example
    /// ```ignore
    /// use algebraic_immunity::vandermonde::VanderMonde;
    ///
    /// let mut mat = VanderMonde::new(vec![vec![1,0,1,0]]);
    /// let vec_to_add = vec![1,0,0,0];
    /// mat.append_row(vec_to_add);
    /// let new_mat = VanderMonde::new(vec![vec![1,0,1,0], vec![1,0,0,0]]);
    /// assert_eq!(mat.elements, new_mat.elements);
    /// ```
    fn append_row(&mut self, v: Vec<u8>) {
        self.elements.push(v)
    }

    fn append_column(&mut self, v: Vec<u8>) {
        for i in 0..self.nrows() {
            self.elements[i].push(v[i]);
        }
    }

    pub fn compute_next(
        &self,
        monom_slice: Vec<usize>,
        support_slice: Vec<usize>,
        idx: usize,
        operations: &Vec<(usize, usize)>,
    ) -> Self {
        let mut m_copy = self.clone();
        let row: Vec<u8> = (0..=idx)
            .map(|i| monomo_eval(support_slice[support_slice.len() - 1], monom_slice[i]) as u8)
            .collect();
        let column: Vec<u8> = (0..idx)
            .map(|i| monomo_eval(support_slice[i], monom_slice[monom_slice.len() - 1]) as u8)
            .collect();

        let n_vect: Vec<u8> = apply_operations(&operations, column);
        m_copy.append_column(n_vect);
        m_copy.append_row(row);

        m_copy
    }

    pub fn compute_vandermonde(support: Vec<usize>, monomials: Vec<usize>) -> Self {
        let result: Vec<Vec<u8>> = support
            .iter()
            .map(|zi| monomials.iter().map(|ej| monomo_eval(*zi, *ej)).collect())
            .collect();
        Self::new(result)
    }

    pub fn fill_rows(&self, support_slice: Vec<usize>, monom_slice: Vec<usize>) -> Self {
        let mut m_copy = self.clone();
        for j in 0..support_slice.len() {
            let row: Vec<u8> = (0..monom_slice.len())
                .map(|i| monomo_eval(support_slice[j], monom_slice[i]) as u8)
                .collect();
            m_copy.append_row(row)
        }

        m_copy
    }

    pub fn construct_and_add_column(
        &self,
        support: &Vec<usize>,
        monom: usize,
        operations: &Vec<(usize, usize)>,
    ) -> Self {
        let mut m_copy = self.clone();
        let column: Vec<u8> = (0..m_copy.nrows())
            .map(|i| monomo_eval(support[i], monom))
            .collect();
        let n_vect: Vec<u8> = apply_operations(&operations, column);
        m_copy.append_column(n_vect);

        m_copy
    }
}

/// Computes the monomial x^u where x and u are elements of F_2^n.
///
/// monomo_eval(x,u) = \sum_{i=0}^n x_i^{u_i}
///
/// # Arguments
///
/// * `s1` - A reference to the element x
/// * `s2` - A reference to the monomial exponent u
///
/// # Returns
///
/// The monomial \sum_{i=0}^n x_i^{u_i}.
///
/// Example
/// ```
/// use algebraic_immunity::vandermonde::monomo_eval;
///
/// assert_eq!(monomo_eval(0b101, 0b010), 0);
/// assert_eq!(monomo_eval(0b011, 0b010), 1);
/// ```
pub fn monomo_eval(x: usize, u: usize) -> u8 {
    ((x & u) == u) as u8
}

fn apply_operations(operations: &Vec<(usize, usize)>, v: Vec<u8>) -> Vec<u8> {
    let mut result = v.clone();
    for &(op1, op2) in operations.iter() {
        result[op1] = (result[op1] + result[op2]) % 2;
    }
    result
}

fn is_submonomial(sub_monom: usize, monom: usize) -> bool {
    (sub_monom & monom) == sub_monom
}

/// Verifies whether a Boolean function `g`, represented in Algebraic Normal Form (ANF),
/// vanishes at all given points in `z`.
///
/// # Arguments
///
/// * `z` - A list of integers, each representing a point `x \in F_2^n`.
/// * `g` - A vector of `0` or `1` coefficients corresponding to the monomials in the ANF.
/// * `mapping` - A vector of integers, each representing the exponent vector `u`
///   of a monomial in the ANF.
///
/// # Returns
///
/// A tuple:
/// * The first element is `true` if `g(x) = 0` for all `x \in z`, otherwise `false`.
/// * The second element is `None` if all evaluations are zero;
///   otherwise, it is `Some((i, x))` where `i` is the index in `z` and `x` is the violating point.
///
/// # Examples
///
/// The function does not vanish in z
/// ```
/// use algebraic_immunity::vandermonde::verify;
///
/// let z = vec![0b110, 0b101];
/// let g = vec![1, 0, 1];  // coefficients for monomials
/// let mapping = vec![0b000, 0b100, 0b110];
/// let result = verify(&z, &g, &mapping);
/// assert_eq!(result.0, false);
/// assert_eq!(result.1.unwrap(), (1, 0b101));
/// ```
///
/// The function vanishes in z
/// ```
/// use algebraic_immunity::vandermonde::verify;
///
/// let z = vec![0b110];
/// let g = vec![1, 0, 1];  // coefficients for monomials
/// let mapping = vec![0b000, 0b100, 0b110];
/// let result = verify(&z, &g, &mapping);
/// assert_eq!(result.0, true);
/// ```
pub fn verify(z: &Vec<usize>, g: &Vec<u8>, mapping: &Vec<usize>) -> (bool, Option<(usize, usize)>) {
    for (idx, &item) in z.iter().enumerate() {
        let anf: Vec<u8> = (0..g.len())
            .filter(|&i| is_submonomial(mapping[i], item))
            .map(|i| g[i])
            .collect();

        if anf.iter().copied().sum::<u8>() % 2 == 1 {
            return (false, Some((idx, item)));
        }
    }
    (true, None)
}
