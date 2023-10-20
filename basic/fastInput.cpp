char buf[1 << 16], *p1 = buf, *p2 = buf;
char get() {
  if (p1 == p2) {
    p1 = buf;
    p2 = p1 + fread(buf, 1, sizeof(buf), stdin);
  }
  if (p1 == p2)
    return -1;
  return *p1++;
}
char readChar() {
  char c = get();
  while (isspace(c))
    c = get();
  return c;
}
int readInt() {
  int x = 0;
  char c = get();
  while (!isdigit(c))
    c = get();
  while (isdigit(c)) {
    x = 10 * x + c - '0';
    c = get();
  }
  return x;
}
