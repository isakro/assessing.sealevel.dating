buf <- sites_sa %>%
  filter(st_is_empty(.) == FALSE) %>%
  st_buffer(dist = 1000)

# hoydedata.no only accepts 10 or less polygons per shape.

tm_shape(st_union(buf, by_feature = FALSE)) +
  tm_polygons()

buffer <- st_union(buf, by_feature = FALSE)
st_write(buf[1,], here::here("analysis/buffer1.shp"))
nrow(buf)


# Function from @TimSalabim copy+pasted directly from:
# https://stackoverflow.com/a/51300037
st_snap_points = function(x, y, max_dist = 60000) {

  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)

  out = do.call(c,
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  return(out)
}
