# Minimal Shapefile fixture for testing.
#
# Creates {dir}/fixture.{shp,shx,dbf,prj} containing three simple CW squares
# with DBF columns NAME (String, width 8), CODE (String, width 2), VAL (Int, width 5).
#
#   idx  NAME   CODE  VAL
#   1    Alpha  A1    1      square at x∈[0,2], y∈[0,2]
#   2    Beta   B2    2      square at x∈[3,5], y∈[0,2]
#   3    Gamma  C3    3      square at x∈[6,8], y∈[0,2]
#
# CRS: EPSG:4326 (WGS 84 geographic).  Coordinates are in degrees but chosen
# small enough that Proj.jl reprojections are also testable.

function _create_test_shapefile(dir::AbstractString) :: String
  mkpath(dir)
  base = joinpath(dir, "fixture")

  # CW squares (ESRI exterior convention: CW = negative signed area).
  rings = [
    [(0.0,0.0),(0.0,2.0),(2.0,2.0),(2.0,0.0),(0.0,0.0)],
    [(3.0,0.0),(3.0,2.0),(5.0,2.0),(5.0,0.0),(3.0,0.0)],
    [(6.0,0.0),(6.0,2.0),(8.0,2.0),(8.0,0.0),(6.0,0.0)],
  ]
  names = ["Alpha", "Beta",  "Gamma"]
  codes = ["A1",    "B2",    "C3"  ]
  vals  = [1,       2,       3     ]

  _fixture_write_shp_shx(base, rings)
  _fixture_write_dbf(base, names, codes, vals)
  write(base * ".prj",
    """GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",""" *
    """SPHEROID["WGS_1984",6378137.0,298.257223563]],""" *
    """PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]""")

  return base * ".shp"
end

# ---------------------------------------------------------------------------
# .shp / .shx
# ---------------------------------------------------------------------------

function _fixture_write_shp_shx(base, rings)
  n = length(rings)

  # Content bytes per record:
  # 4 (type) + 32 (bbox) + 4 (num_parts) + 4 (num_pts) + 4 (parts[0]) + npts*16
  cb = [4 + 32 + 4 + 4 + 4 + length(r)*16 for r in rings]
  cw = cb .÷ 2   # content in 16-bit words

  shp_words = 50 + sum(4 + w for w in cw)   # 50 = 100-byte header
  shx_words = 50 + 4*n

  all_x = [p[1] for r in rings for p in r]
  all_y = [p[2] for r in rings for p in r]
  bbox  = (minimum(all_x), minimum(all_y), maximum(all_x), maximum(all_y))

  open(base * ".shp", "w") do shp
    open(base * ".shx", "w") do shx
      _fixture_shp_header(shp, shp_words, bbox)
      _fixture_shp_header(shx, shx_words, bbox)

      offset = 50   # running shp offset in 16-bit words
      for (i, ring) in enumerate(rings)
        # SHX entry.
        write(shx, hton(Int32(offset)))
        write(shx, hton(Int32(cw[i])))

        # SHP record header (big-endian).
        write(shp, hton(Int32(i)))
        write(shp, hton(Int32(cw[i])))

        # SHP record content (little-endian).
        npts = length(ring)
        rx   = [p[1] for p in ring]
        ry   = [p[2] for p in ring]
        write(shp, htol(Int32(5)))                   # shape type: Polygon
        write(shp, htol(Float64(minimum(rx))))
        write(shp, htol(Float64(minimum(ry))))
        write(shp, htol(Float64(maximum(rx))))
        write(shp, htol(Float64(maximum(ry))))
        write(shp, htol(Int32(1)))                   # num_parts
        write(shp, htol(Int32(npts)))                # num_points
        write(shp, htol(Int32(0)))                   # parts[0]
        for (x, y) in ring
          write(shp, htol(Float64(x)))
          write(shp, htol(Float64(y)))
        end

        offset += 4 + cw[i]
      end
    end
  end
end

function _fixture_shp_header(io, length_words, bbox)
  write(io, hton(Int32(9994)))           # file code
  for _ in 1:5; write(io, hton(Int32(0))); end
  write(io, hton(Int32(length_words)))   # file length in 16-bit words
  write(io, htol(Int32(1000)))           # version
  write(io, htol(Int32(5)))             # shape type: Polygon
  xmin, ymin, xmax, ymax = bbox
  for v in (xmin, ymin, xmax, ymax, 0.0, 0.0, 0.0, 0.0)
    write(io, htol(Float64(v)))
  end
end

# ---------------------------------------------------------------------------
# .dbf (dBASE III+)
# ---------------------------------------------------------------------------

function _fixture_write_dbf(base, names, codes, vals)
  n        = length(names)
  fw_name  = 8
  fw_code  = 2
  fw_val   = 5
  rec_size = 1 + fw_name + fw_code + fw_val
  hdr_size = 32 + 3*32 + 1   # file hdr + 3 field descriptors + 0x0D terminator

  open(base * ".dbf", "w") do io
    # File header (32 bytes).
    write(io, UInt8(0x03))            # dBASE III+ version
    write(io, UInt8(126))             # year since 1900 (2026)
    write(io, UInt8(3))               # month
    write(io, UInt8(16))              # day
    write(io, htol(Int32(n)))
    write(io, htol(Int16(hdr_size)))
    write(io, htol(Int16(rec_size)))
    write(io, zeros(UInt8, 20))

    # Field descriptors.
    _fixture_dbf_field(io, "NAME", 'C', fw_name)
    _fixture_dbf_field(io, "CODE", 'C', fw_code)
    _fixture_dbf_field(io, "VAL",  'N', fw_val)
    write(io, UInt8(0x0D))

    # Records.
    for i in 1:n
      write(io, UInt8(' '))
      write(io, rpad(names[i], fw_name)[1:fw_name])
      write(io, rpad(codes[i], fw_code)[1:fw_code])
      write(io, lpad(string(vals[i]), fw_val)[1:fw_val])
    end
    write(io, UInt8(0x1A))   # EOF
  end
end

function _fixture_dbf_field(io, name, type, width)
  nb = zeros(UInt8, 11)
  for (i, c) in enumerate(name); i > 11 && break; nb[i] = UInt8(c); end
  write(io, nb)
  write(io, UInt8(type))
  write(io, zeros(UInt8, 4))
  write(io, UInt8(width))
  write(io, UInt8(0))
  write(io, zeros(UInt8, 14))
end
