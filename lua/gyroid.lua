function surface (x, y, z)
  return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x) - params.isolevel
end