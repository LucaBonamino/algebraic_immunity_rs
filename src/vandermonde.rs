use lin_algebra::gf2_matrix::GF2Matrix;
use lin_algebra::matrix::MatrixTrait;
use std::ops::{Deref, DerefMut};

#[derive(Clone)]
pub struct VanderMonde{
    pub matrix: GF2Matrix
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


impl VanderMonde{

    pub fn new(elements: Vec<Vec<u8>>) -> Self{
        VanderMonde { matrix: GF2Matrix::new(elements) }
    }

    pub fn __repr__(&self) -> String {
        let rows: Vec<String> = self
            .elements
            .iter()
            .map(|row| format!("{:?}", row))
            .collect();
        format!("[{}]", rows.join(", "))
    }


    fn append_row(&mut self, v: Vec<u8>) {
        self.elements.push(v)
    }

    fn append_column(&mut self, v: Vec<u8>) {
        for i in 0..self.nrows() {
            self.elements[i].push(v[i]);
        }
    }

    pub fn rank(&self) -> usize {
        let mut count = 0;
        let mut pivot_columns = std::collections::HashSet::new();

        for i in 0..self.nrows() {
            let p = GF2Matrix::get_pivot(&self.elements[i]);
            if let Some(col) = p {
                if pivot_columns.insert(col) {
                    count += 1;
                }
            }
        }
        count
    }


    pub fn compute_next(
        &self,
        monom_slice: Vec<String>,
        support_slice: Vec<String>,
        idx: usize,
        operations: Vec<(usize, usize)>
    ) -> Self {
        let mut m_copy = self.clone();
        let row: Vec<u8> = (0..=idx)
            .map(|i| str_ops(&support_slice[support_slice.len() - 1], &monom_slice[i]) as u8)
            .collect();
        let column: Vec<u8> = (0..idx)
            .map(|i| str_ops(&support_slice[i], &monom_slice[monom_slice.len() - 1]) as u8)
            .collect();

        let n_vect: Vec<u8> = apply_operations(&operations, column);
        m_copy.append_column(n_vect);
        m_copy.append_row(row);

        m_copy
    }

}


pub fn str_ops(s1: &str, s2: &str) -> u8 {
    s1.chars()
        .zip(s2.chars())
        .map(|(c1, c2)| {
            let base = c1.to_digit(10).unwrap() as u8;
            let exp = c2.to_digit(10).unwrap() as u8;
            base.pow(exp as u32)
        })
        .product()
}

fn apply_operations(operations: &Vec<(usize, usize)>, v: Vec<u8>) -> Vec<u8> {
    let mut result = v.clone();
    for &(op1, op2) in operations.iter() {
        result[op1] = (result[op1] + result[op2]) % 2;
    }
    result
}


fn is_submonomial(sub_monom: &str, monom: &str) -> bool {
    assert_eq!(sub_monom.len(), monom.len(), "The lengths of sub_monom and monom must be equal");

    for (char1, char2) in sub_monom.chars().zip(monom.chars()) {
        if char1 > char2 {
            return false;
        }
    }
    true
}

pub fn verify(z: Vec<String>, g: Vec<u8>, mapping: Vec<String>) -> (bool, Option<(usize, String)>) {
    for (idx, item) in z.iter().enumerate() {
        let anf: Vec<u8> = (0..g.len())
            .filter(|&i| is_submonomial(&mapping[i], item))
            .map(|i| g[i])
            .collect();

        if anf.iter().copied().sum::<u8>() % 2 == 1 {
            return (false, Some((idx, item.clone())));
        }
    }
    (true, None)
}


