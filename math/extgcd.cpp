array<i64, 3> extgcd(i64 a, i64 b) {
	if (!b) { return {a, 1, 0}; }
	auto [g, x, y] = extgcd(b, a % b);
	return {g, y, x - a / b * y};
}
