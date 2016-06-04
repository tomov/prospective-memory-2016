function d = CohensD(m1, s1, n1, m2, s2, n2)

s = sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
d = (m1 - m2) / s;

end

