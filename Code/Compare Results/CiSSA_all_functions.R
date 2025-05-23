circulant_SSA <- function (ts, L = NULL, extend_flag = FALSE) 
{
    time_series <- ts
    N <- length(time_series)
    if (is.null(L)) {
        L <- (N + 1)%/%2
    }
    if (extend_flag == FALSE) {
        H <- 0
        time_series <- ts
    }
    else {
        H <- L
        time_series <- extend(ts, H)
    }
    X <- hankel(time_series, L)
    if (L%%2) {
        nf2 <- (L + 1)/2 - 1
    }
    else {
        nf2 <- L/2 - 1
    }
    nft <- nf2 + abs((L%%2) - 2)
    autocov <- numeric(L)
    for (m in 0:(L - 1)) {
        autocov[[m + 1]] <- sum(time_series[1:(N - m)] * time_series[(1 + 
            m):N])/(N - m)
    }
    circ_first_row <- numeric(L)
    for (m in 0:(L - 1)) {
        circ_first_row[[m + 1]] <- (L - m)/L * autocov[[m + 1]] + 
            (m)/L * autocov[[L - m]]
    }
    S_C <- matrix(circ_first_row, nrow = 1)
    shifted_vector <- circ_first_row
    for (i in 2:(L)) {
        shifted_vector <- shift_vector(shifted_vector)
        S_C <- rbind(as.vector(shifted_vector), S_C)
    }
    U <- dftmtx(L)/sqrt(L)
    U[, 1] <- Re(U[, 1])
    for (k in 1:nf2) {
        u_k <- U[, k + 1]
        U[, k + 1] <- sqrt(2) * Re(u_k)
        U[, L + 2 - (k + 1)] <- sqrt(2) * Im(u_k)
    }
    if (L%%2 != 0) {
        U[, nft] <- Re(U[, nft])
    }
    psd <- abs(diag(t(U) %*% S_C %*% U))
    W <- t(U) %*% X
    R <- matrix(0, nrow = N + 2 * H, ncol = L)
    for (k in 1:L) {
        R[, k] <- diag_averaging(U[, k] %*% t(W[k, ]))
    }
    Z <- matrix(0, nrow = N + 2 * H, ncol = nft)
    Z[, 1] <- R[, 1]
    imp <- numeric(nft)
    lambda_sm <- sum(psd)
    imp[1] <- psd[1]/lambda_sm
    for (k in 1:nf2) {
        Z[, k + 1] <- R[, k + 1] + R[, L + 2 - (k + 1)]
        imp[k + 1] <- (psd[k + 1] + psd[L + 2 - (k + 1)])/lambda_sm
    }
    if (L%%2 != 0) {
        Z[, nft] <- R[, nft]
        imp[nft] <- psd[nft]/lambda_sm
    }
    list(t_series = Z[(H + 1):(N + H), ], importance = imp, freq = (0:(length(imp) - 
        1))/L)
}


cissa_like_fourier_transform <- function (ts, L) 
{
    reconstruct_fft <- function(x, y, frequencies) {
        fft_y <- fft(y)
        amplitudes <- Mod(fft_y)
        phases <- Arg(fft_y)
        reconstructed <- matrix(0, length(amplitudes), length(x))
        n <- length(amplitudes)
        for (i in 1:(length(amplitudes))) {
            reconstructed[i, ] <- amplitudes[i] * cos(2 * pi * 
                frequencies[i] * (x) + phases[i])/n
        }
        return(list(ts = reconstructed, frequencies = frequencies))
    }
    n <- length(ts)
    x <- 0:(n - 1)
    L <- L
    frequencies <- (0:(L - 1))/L
    y <- ts
    y_main <- y
    x_main <- x
    X <- hankel(y, L)
    K <- dim(X)[2]
    res <- list()
    for (i in 1:K) {
        y <- X[, i]
        x <- x_main[1:(1 + L - 1)]
        res[[i]] <- reconstruct_fft(x, y, frequencies)$ts
    }
    res_mult <- res
    averaging <- function(res_comp_wise_mult) {
        K <- dim(X)[2]
        counters <- rep(0, n)
        res <- matrix(0, ncol = n, nrow = L)
        for (i in 1:length(res_comp_wise_mult)) {
            res[, i:(i + L - 1)] <- res[, i:(i + L - 1)] + res_comp_wise_mult[[i]]
            counters[i:(i + L - 1)] <- counters[i:(i + L - 1)] + 
                1
        }
        for (i in 1:n) {
            res[, i] <- res[, i]/counters[i]
        }
        res
    }
    avr <- averaging(res_mult)
    group_by_elementary_freq_foureir <- function(res_averaged) {
        nf2 <- 0
        if (L%%2) {
            nf2 <- (L + 1)/2 - 1
        }
        else {
            nf2 <- L/2 - 1
        }
        nft <- nf2 + abs((L%%2) - 2)
        Z <- matrix(0, ncol = nft, nrow = n)
        Z[, 1] <- res_averaged[1, ]
        for (k in 1:nf2) {
            Z[, k + 1] <- res_averaged[k + 1, ] + res_averaged[L + 
                2 - (k + 1), ]
        }
        if (L%%2 != 0) {
            Z[, nft] <- res_averaged[nft, ]
        }
        return(list(t_series = Z, freq = (0:dim(Z)[2])/L))
    }
    rs <- group_by_elementary_freq_foureir(avr)
    return(rs)
}


clust_basis_general <- function (U, roots, type = c("hierarchial", "distance", "frequency"), 
    k = 2, auto_trend_freq = 1/24, delta = 0.001, h = NULL, order = FALSE) 
{
    ord <- order(abs(Arg(roots)))
    roots <- roots[ord]
    U <- U[, ord, drop = FALSE]
    stopifnot(length(k) == 1)
    maxk <- sum(Im(roots) >= -.Machine$double.eps)
    if ((length(roots) == 1) | (k > maxk)) 
        k <- 1
    groups <- numeric(0)
    if (type == "hierarchial") {
        d <- stats::dist(cbind(Re(roots), abs(Im(roots))), method = "euclidian")
        hc <- hclust(d, method = "complete")
        idx <- cutree(hc, k = k, h = h)
        groups <- tapply(seq_along(idx), idx, identity)
        names(groups) <- paste("F", names(groups), sep = "")
    }
    else if (type == "distance") {
        groups <- dist_clust(roots, delta = delta)
    }
    else {
        groups <- freq_clust(roots, auto_trend_freq)
    }
    for (group in groups) {
        group_vecs <- U[, group, drop = FALSE]
        U[, group] <- svd(cbind(Re(group_vecs), Im(group_vecs)), 
            nu = length(group), nv = 0)$u
    }
    U <- Re(U)
    if (order) {
        U <- U[, unlist(groups), drop = FALSE]
        l <- sapply(groups, length)
        cl <- cumsum(c(0, l))
        groups <- lapply(seq_along(l), function(i) (cl[i] + 1):cl[i + 
            1])
    }
    list(basis = U, groups = groups)
}


df_to_latex_print <- function (df) 
{
    df[] <- lapply(df, function(x) {
        if (is.numeric(x)) {
            formatC(x, format = "e", digits = 1)
        }
        else {
            x
        }
    })
    table_latex <- xtable(df, caption = "Example Table")
    print(table_latex, include.rownames = FALSE)
}


dftmtx <- function (n) 
{
    y <- stats::mvfft(diag(1, n))
    y
}


diag_averaging <- function (A) 
{
    B <- Re(A[nrow(A):1, ])
    as.numeric(lapply(split(B, -(row(B) - col(B))), mean))
}


dist_clust <- function (roots, delta) 
{
    roots_data <- cbind(Re(roots), abs(Im(roots)))
    split_mask <- (roots_data[, 1] < 0) & (roots_data[, 2] < 
        0.1) & roots_data[, 2] > -0.1
    part1 <- roots_data[split_mask, ]
    if (!is.matrix(part1)) 
        part1 <- matrix(part1, nrow = 1, ncol = 2)
    part2 <- roots_data[split_mask == FALSE, ]
    if (!is.matrix(part2)) 
        part2 <- matrix(part2, nrow = 1, ncol = 2)
    if (length(part1) > 0) 
        indices1 <- dist_clust_table(part1, delta = delta)
    else indices1 <- numeric(0)
    indices2 <- dist_clust_table(part2, delta = delta)
    if (length(indices1) > 0) 
        max_index1 <- max(indices1)
    else max_index1 <- 0
    indices2 <- indices2 + max_index1
    res_indices <- numeric(length(roots))
    res_indices[split_mask] <- indices1
    res_indices[split_mask == FALSE] <- indices2
    groups <- list()
    res_k <- max(res_indices)
    for (i in 1:res_k) {
        clust_indices <- which(res_indices == i)
        groups[[paste0("F", i)]] <- clust_indices
    }
    return(groups)
}


dist_clust_old <- function (roots, delta) 
{
    pdf(file = "tex/fig/signal_roots.pdf", width = 6.5, height = 4)
    plot(roots, asp = 1)
    dev.off()
    roots_data <- cbind(Re(roots), abs(Im(roots)))
    res_indices <- numeric(0)
    res_k <- numeric(0)
    max_k <- nrow(unique(roots_data))
    if (max_k == 1) {
        return(list(F1 = 1:nrow(roots_data)))
    }
    k <- 2
    while (k <= min(max_k, nrow(roots_data) - 1)) {
        res <- kmeans(roots_data, centers = k, nstart = 10)
        if (res$tot.withinss/res$totss < delta) {
            res_k <- k
            res_indices <- res$cluster
            break
        }
        k <- k + 1
    }
    if (length(res_indices) == 0) {
        res_k <- nrow(roots_data)
        res_indices <- 1:res_k
    }
    groups <- list()
    for (i in 1:res_k) {
        clust_indices <- which(res_indices == i)
        groups[[paste0("F", i)]] <- clust_indices
    }
    return(groups)
}


dist_clust_table <- function (roots_table, delta) 
{
    set.seed(10)
    roots_data <- roots_table
    res_indices <- numeric(0)
    res_k <- numeric(0)
    if (is.null(ncol(roots_table))) 
        roots_data <- matrix(roots_data, nrow = 1)
    max_k <- nrow(unique(roots_data))
    if (max_k == 0) {
        return(numeric(0))
    }
    if (max_k == 1) {
        return(rep(1, nrow(roots_table)))
    }
    k <- 2
    while (k <= min(max_k, nrow(roots_data) - 1)) {
        res <- kmeans(roots_data, centers = k, nstart = 10)
        if (res$tot.withinss/res$totss < delta) {
            res_k <- k
            res_indices <- res$cluster
            break
        }
        k <- k + 1
    }
    if (length(res_indices) == 0) {
        res_k <- nrow(roots_data)
        if (res_k == 1) 
            res_indices <- 1
        else res_indices <- 1:res_k
    }
    return(res_indices)
}


do_one_experiment <- function (methods_list, methods_names) 
{
}


eossa_new <- function (x, nested.groups, clust_type = c("hierarchial", "distance", 
    "frequency"), k = 2, delta = 0.001, auto_trend_freq = 1/24, 
    subspace = c("column", "row"), dimensions = NULL, solve.method = c("ls", 
        "tls"), beta = 8, ...) 
{
    if (missing(nested.groups)) 
        nested.groups <- as.list(1:min(nsigma(x), nu(x)))
    subspace <- match.arg(subspace)
    solve.method <- match.arg(solve.method)
    Rssa:::.maybe.continue(x, groups = nested.groups, ...)
    idx <- sort(unique(unlist(nested.groups)))
    triples <- Rssa:::.get.orth.triples(x, idx, do.orthogonalize = FALSE)
    osigma <- triples$sigma
    U <- triples$U
    V <- triples$V
    if (identical(subspace, "column")) {
        vectors <- U
        V <- V * rep(osigma, each = nrow(V))
        wmask <- Rssa:::.wmask(x)
    }
    else if (identical(subspace, "row")) {
        vectors <- V
        U <- U * rep(sqrt(osigma), each = nrow(U))
        wmask <- Rssa:::.fmask(x)
    }
    sigma <- rep(1, length(idx))
    if (is.null(dimensions)) {
        dimensions <- seq_len(Rssa:::.dim(x))
    }
    d <- dim(wmask)
    if (max(dimensions) > length(d)) {
        stop(sprintf("some of input dimension indices exceed the actual number of object dimensions (%d)", 
            length(d)))
    }
    Zs <- lapply(dimensions, function(ndim) {
        Rssa:::.shift.matrix(vectors, wmask = wmask, ndim = ndim, 
            circular = x$circular[ndim], solve.method = solve.method)
    })
    Z <- Rssa:::.matrix.linear.combination(Zs, beta)
    Ze <- eigen(Z, symmetric = FALSE)
    mb <- clust_basis_general(Ze$vectors, Ze$values, type = clust_type, 
        k = k, delta = delta, auto_trend_freq = auto_trend_freq)
    C <- mb$basis
    nested.groups <- mb$groups
    U <- U %*% C
    V <- t(qr.solve(C, t(V)))
    x <- clone(x, copy.cache = FALSE)
    Rssa:::.save.oblique.decomposition(x, sigma, U, V, idx)
    nested.groups <- lapply(nested.groups, function(group) idx[group])
    iossa.groups.all <- Rssa:::.get(x, "iossa.groups.all", allow.null = TRUE)
    if (is.null(iossa.groups.all)) {
        iossa.groups.all <- list()
    }
    valid.groups <- as.logical(sapply(iossa.groups.all, function(group) length(intersect(group, 
        idx)) == 0))
    Rssa:::.set(x, "iossa.groups", nested.groups)
    Rssa:::.set(x, "iossa.groups.all", c(nested.groups, iossa.groups.all[valid.groups]))
    Rssa:::.set(x, "ossa.set", idx)
    if (!is.null(Rssa:::.decomposition(x, "nPR"))) {
        if (any(idx <= Rssa:::.decomposition(x, "nPR"))) {
            Rssa:::.set.decomposition(x, nPR = 0, nPL = 0)
        }
        else if (any(idx <= sum(unlist(Rssa:::.decomposition(x, 
            c("nPR", "nPL")))))) {
            Rssa:::.set.decomposition(x, nPL = 0)
        }
    }
    if (!inherits(x, "ossa")) {
        class(x) <- c("ossa", class(x))
    }
    x$call <- match.call()
    invisible(x)
}


extend <- function (x, H) 
{
    N <- length(x)
    p <- floor(N/3)
    dx <- diff(x)
    A <- aryule(dx, p)$a
    y <- x
    dy <- diff(y)
    er <- signal::filter(A, 1, dy)
    dy <- signal::filter(1, A, c(er, rep(0, H)))
    y <- y[1] + c(0, cumsum(dy))
    y <- rev(y)
    dy <- diff(y)
    er <- signal::filter(A, 1, dy)
    dy <- signal::filter(1, A, c(er, rep(0, H)))
    y <- y[1] + c(0, cumsum(dy))
    xe <- rev(y)
    xe
}


f_const <- function (x, C = 0) 
{
    rep(C, length(x))
}


f_cos <- function (x, A = 1, omega = 1/4, phi = 0) 
{
    f_exp_mod_harm_series(x, A, alpha = 0, omega = omega, phi = phi)
}


f_exp <- function (x, A = 1, alpha = 1) 
{
    A * exp(alpha * x)
}


f_exp_cos <- function (x, A = 1, alpha = 1, omega = 1/4, phi = 0) 
{
    f_exp_mod_harm_series(x, A, alpha, omega, phi)
}


f_exp_mod_harm_series <- function (x, A = 1, alpha = 1, omega = 1/4, phi = 0) 
{
    A * exp(alpha * x) * cos(2 * pi * omega * x + phi)
}


f_linear <- function (x, a = 1, b = 0) 
{
    a * x + b
}


f_sin <- function (x, A = 1, omega = 1/4, phi = 3 * pi/2) 
{
    f_exp_mod_harm_series(x, A, alpha = 0, omega = omega, phi = phi)
}


f_sum <- function (x) 
{
    f_cos(x, omega = omega_cs) + f_exp_mod_harm_series(x, a = a, 
        omega = omega_exp) + f_sin(x, omega = omega_sn)
}


freq_clust <- function (roots, auto_trend_freq = 1/24) 
{
    freqs <- abs(Im(log(roots)))/2/pi
    trend_group <- which(freqs < auto_trend_freq)
    periodic_group <- which(freqs >= auto_trend_freq)
    if (length(trend_group) == 0) 
        return(list(F1 = periodic_group))
    if (length(periodic_group) == 0) 
        return(list(F1 = trend_group))
    return(list(F1 = trend_group, F2 = periodic_group))
}


generate_ts <- function (func, n = 1000, ...) 
{
    ts(func(1:n, ...))
}


good_experiment <- function (L, n, x, t_series, groups_ssa, eossa_nested_groups, 
    groups_freqs, table, y = NULL) 
{
    if (is.null(y)) {
        y <- Reduce("+", t_series)
    }
    y_real <- Reduce("+", t_series)
    s_ssa <- ssa(y[1:(n - 1)], L)
    r_ssa <- reconstruct(s_ssa, groups = groups_ssa)
    e_ssa <- eossa_new(s_ssa, nested.groups = eossa_nested_groups, 
        clust_type = "distance")
    g_sesonal_e <- grouping.auto(e_ssa, base = "eigen", freq.bins = groups_freqs, 
        threshold = 0.5)
    r_ssa_e <- reconstruct(e_ssa, groups = g_sesonal_e)
    r_fft <- reconstruct_fft(x, y)
    r_fft_grouped <- grouping_cissa(r_fft, groups = groups_freqs)$t_series
    r_fft_extended <- reconstruct_fft(x, y, TRUE)
    r_fft_grouped_extended <- grouping_cissa(r_fft_extended, 
        groups = groups_freqs)$t_series
    r_cissa <- circulant_SSA(y, L)
    r_cissa_grouped <- grouping_cissa(r_cissa, groups = groups_freqs)$t_series
    r_cissa_ext <- circulant_SSA(y, L, extend_flag = TRUE)
    r_cissa_grouped_ext <- grouping_cissa(r_cissa_ext, groups = groups_freqs)$t_series
    ssa_sum <- 0
    ssa_e_sum <- 0
    fft_sum <- 0
    fft_ext_sum <- 0
    cissa_sum <- 0
    cissa_ext_sum <- 0
    for (i in 1:length(t_series)) {
        ts <- t_series[[i]]
        table[1, i + 1] <- mse(ts[1:(n - 1)], r_ssa[[i]])
        ssa_sum <- ssa_sum + r_ssa[[i]]
        table[2, i + 1] <- mse(ts[1:(n - 1)], r_ssa_e[[i]])
        ssa_e_sum <- ssa_e_sum + r_ssa_e[[i]]
        table[3, i + 1] <- mse(ts, r_fft_grouped[[i]])
        fft_sum <- fft_sum + r_fft_grouped[[i]]
        table[4, i + 1] <- mse(ts, r_fft_grouped_extended[[i]])
        fft_ext_sum <- fft_ext_sum + r_fft_grouped_extended[[i]]
        table[5, i + 1] <- mse(ts, r_cissa_grouped[[i]])
        cissa_sum <- cissa_sum + r_cissa_grouped[[i]]
        table[6, i + 1] <- mse(ts, r_cissa_grouped_ext[[i]])
        cissa_ext_sum <- cissa_ext_sum + r_cissa_grouped_ext[[i]]
    }
    table[1, length(t_series) + 2] <- mse(y_real[1:(n - 1)], 
        ssa_sum)
    table[2, length(t_series) + 2] <- mse(y_real[1:(n - 1)], 
        ssa_e_sum)
    table[3, length(t_series) + 2] <- mse(y_real, fft_sum)
    table[4, length(t_series) + 2] <- mse(y_real, fft_ext_sum)
    table[5, length(t_series) + 2] <- mse(y_real, cissa_sum)
    table[6, length(t_series) + 2] <- mse(y_real, cissa_ext_sum)
    return(table)
}



grouping_cissa <- function (cissa_res, groups) 
{
  freq <- cissa_res$freq
  t_series <- cissa_res$t_series
  residuals <- 0
  result <- setNames(as.list(rep(0, length(groups))), names(groups))
  result_freqs <- list()
  for (i in 1:length(cissa_res$freq)) {
    flag <- FALSE
    for (name in names(groups)) {
      if (flag == TRUE) {
        break
      }
      if (groups[[name]][1] <= freq[i] & freq[i] <= groups[[name]][2]) {
        flag <- TRUE
        result[[name]] <- result[[name]] + t_series[, 
                                                    i]
        result_freqs[[name]] <- c(result_freqs[[name]], 
                                  freq[i])
      }
    }
    if (flag == FALSE) {
      residuals <- residuals + t_series[, i]
    }
  }
  result[["residuals"]] <- residuals
  return(list(t_series = result, freqs_by_group = result_freqs))
  result
}



interesting_sin <- function (x) 
{
    first_seq <- sin(x[x <= 15] * 2 * pi/12)
    plot(sin(x[x <= 15] * 2 * pi/12), type = "l")
    second_seq <- x[15 < x & x <= 45]/15 * sin(x[15 < x & x <= 
        45] * 2 * pi/12)
    plot(second_seq, type = "l")
    third_seq <- 3 * sin(x[x > 45] * 2 * pi/12)
    c(first_seq, second_seq, third_seq)
}


mse <- function (f_true, f_reconstructed) 
{
    mean((f_true - f_reconstructed)^2)
}


n_mse_tests <- function (n) 
{
    n <- 96 * 2 - 1
    L <- 96
    sigma <- 0.1
    C <- 1
    omega_cs <- 1/12
    omega_sn <- 1/24
    a <- 1/100
    f_sum <- function(x) {
        f_const(x, C = C) + f_cos(x, omega = omega_cs) + f_exp(x, 
            a = a) + f_sin(x, omega = omega_sn)
    }
    f_C <- generate_ts(f_const, n, C = C)
    f_c <- generate_ts(f_cos, n, omega = omega_cs)
    f_s <- generate_ts(f_sin, n, omega = omega_sn)
    f_e <- generate_ts(f_exp, n, a = a)
    mse_lst <- list()
    for (i in 1:n) {
        f_noise <- rnorm(n, sd = sigma)
        f_n <- f_sum(1:n) + f_noise
        c <- circulant_SSA(f_n, L = L, extend_flag = TRUE)
        r <- grouping_cissa(c, groups = list(trend = c(0, 1/1000), 
            sesonal2 = c(1/25, 1/23), sesonal1 = c(1/13, 1/10)))$t_series
        mse_lst$cissa <- c(mse_lst$cissa, mse(f_sum(1:n), r$trend + 
            r$sesonal1 + r$sesonal2))
        s <- ssa(f_n, L)
        e <- fossa(s)
        g_sesonal <- grouping.auto(e, base = "eigen", freq.bins = list(trend = 1/1000, 
            sesonal2 = c(1/25, 1/23), sesonal1 = c(1/13, 1/10)), 
            threshold = 0.5)
        r <- reconstruct(e, groups = c(list(exp = 1, C = 2), 
            g_sesonal))
        mse_lst$ssa <- c(mse_lst$ssa, mse(f_sum(1:n), r$trend + 
            r$sesonal2 + r$sesonal1))
    }
    return(mse_lst)
}


reconstruct_fft <- function (x_init, y_init, extend_flag = FALSE) 
{
    x <- x_init
    y <- y_init
    N <- length(y_init)
    H <- 0
    if (extend_flag == TRUE) {
        H <- length(y)%/%2
        y <- extend(y, H)
        x <- 0:(length(y) - 1)
    }
    frequencies <- (0:(length(x) - 1))/length(x)
    fft_y <- fft(y)
    amplitudes <- Mod(fft_y)
    phases <- Arg(fft_y)
    reconstructed <- matrix(0, length(amplitudes), length(x))
    n <- length(amplitudes)
    L <- n
    for (i in 1:(length(amplitudes))) {
        reconstructed[i, ] <- amplitudes[i] * cos(2 * pi * frequencies[i] * 
            (x) + phases[i])/n
    }
    group_by_elementary_freq_foureir <- function(res_averaged) {
        nf2 <- 0
        if (L%%2) {
            nf2 <- (L + 1)/2 - 1
        }
        else {
            nf2 <- L/2 - 1
        }
        nft <- nf2 + abs((L%%2) - 2)
        Z <- matrix(0, ncol = nft, nrow = n)
        Z[, 1] <- res_averaged[1, ]
        for (k in 1:nf2) {
            Z[, k + 1] <- res_averaged[k + 1, ] + res_averaged[L + 
                2 - (k + 1), ]
        }
        if (L%%2 != 0) {
            Z[, nft] <- res_averaged[nft, ]
        }
        return(list(t_series = Z[(H + 1):(N + H), ], freq = (0:(dim(Z)[2] - 
            1))/L))
    }
    rs <- group_by_elementary_freq_foureir(reconstructed)
    return(list(t_series = rs$t_series, freq = rs$freq))
}


shift_vector <- function (vec) 
{
    last_element <- tail(vec, 1)
    vec <- vec[-length(vec)]
    shifted_vec <- c(last_element, vec)
    return(shifted_vec)
}


