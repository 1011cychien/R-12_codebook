template <typename T>
T minRotation(T s) {
  int n = s.size();
  int i = 0, j = 1;
  s.insert(s.end(), s.begin(), s.end());
  while (i < n && j < n) {
    int k = 0;
    while (k < n && s[i + k] == s[j + k]) {
      k++;
    }
    if (s[i + k] <= s[j + k]) {
      j += k + 1;
    } else {
      i += k + 1;
    }
    if (i == j) {
      j++;
    }
  }
  int ans = i < n ? i : j;
  return T(s.begin() + ans, s.begin() + ans + n);
}
