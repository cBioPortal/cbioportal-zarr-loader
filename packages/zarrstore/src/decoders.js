import * as zarr from "zarrita";

// Default opener tries v2 first, then v3
async function defaultOpen(location, opts) {
  try {
    return await zarr.open.v2(location, opts);
  } catch {
    return await zarr.open.v3(location, opts);
  }
}

export async function readArray(arr) {
  const chunk = await zarr.get(arr);
  return { data: chunk.data, shape: chunk.shape };
}

export function toStringArray(data) {
  if (Array.isArray(data)) return data;
  if (typeof data[0] === "string") return Array.from(data);
  // zarrita may return a TypedArray of bytes for vlen-utf8; decode if needed
  if (data instanceof Uint8Array) {
    const decoder = new TextDecoder();
    return decoder.decode(data).split("\0").filter(Boolean);
  }
  return Array.from(data, (v) => String(v));
}

export async function decodeCategorical(group, open = defaultOpen) {
  const codes = await open(group.resolve("codes"), { kind: "array" });
  const categories = await open(group.resolve("categories"), {
    kind: "array",
  });
  const codesResult = await readArray(codes);
  const categoriesResult = await readArray(categories);
  const ordered = group.attrs?.ordered ?? false;

  const catValues =
    typeof categoriesResult.data[0] === "string" ||
    categoriesResult.data instanceof Array
      ? toStringArray(categoriesResult.data)
      : categoriesResult.data;

  const values = Array.from(codesResult.data, (code) =>
    code < 0 ? null : catValues[code],
  );

  return { values, categories: catValues, ordered };
}

export async function decodeColumn(group, colName, open = defaultOpen) {
  let node;
  try {
    node = await open(group.resolve(colName), { kind: "group" });
  } catch {
    // not a group — open as array
    const arr = await open(group.resolve(colName), { kind: "array" });
    const result = await readArray(arr);
    if (typeof result.data[0] === "string" || result.data instanceof Array) {
      return toStringArray(result.data);
    }
    return result.data;
  }

  const encodingType = node.attrs?.["encoding-type"];
  if (encodingType === "categorical") {
    const decoded = await decodeCategorical(node, open);
    return decoded.values;
  }
  if (
    encodingType === "nullable-integer" ||
    encodingType === "nullable-boolean"
  ) {
    const decoded = await decodeNullable(node, open);
    return decoded.values;
  }
  // fallback: try reading as categorical (common even without explicit encoding-type)
  try {
    const decoded = await decodeCategorical(node, open);
    return decoded.values;
  } catch {
    throw new Error(
      `Unknown column encoding-type "${encodingType}" for column "${colName}"`,
    );
  }
}

export async function decodeDataframe(group, open = defaultOpen) {
  const attrs = group.attrs;
  const indexKey = attrs["_index"];
  const columnOrder = attrs["column-order"];

  const indexArr = await open(group.resolve(indexKey), { kind: "array" });
  const indexResult = await readArray(indexArr);
  const index = toStringArray(indexResult.data);

  const columns = {};
  for (const colName of columnOrder) {
    columns[colName] = await decodeColumn(group, colName, open);
  }

  return { index, columns, columnOrder: Array.from(columnOrder) };
}

export async function decodeNullable(group, open = defaultOpen) {
  const valuesArr = await open(group.resolve("values"), {
    kind: "array",
  });
  const maskArr = await open(group.resolve("mask"), { kind: "array" });
  const valuesResult = await readArray(valuesArr);
  const maskResult = await readArray(maskArr);

  const values = Array.from(valuesResult.data, (v, i) =>
    maskResult.data[i] ? null : v,
  );

  return { values, mask: maskResult.data };
}

export async function decodeSparseMatrix(group, open = defaultOpen) {
  const attrs = group.attrs;
  const format = attrs["encoding-type"]?.replace("_matrix", ""); // "csr" or "csc"
  const shape = attrs.shape;

  const dataArr = await open(group.resolve("data"), { kind: "array" });
  const indicesArr = await open(group.resolve("indices"), {
    kind: "array",
  });
  const indptrArr = await open(group.resolve("indptr"), {
    kind: "array",
  });

  const [dataResult, indicesResult, indptrResult] = await Promise.all([
    readArray(dataArr),
    readArray(indicesArr),
    readArray(indptrArr),
  ]);

  return {
    format,
    data: dataResult.data,
    indices: indicesResult.data,
    indptr: indptrResult.data,
    shape,
  };
}

export function sparseToDense(sparse) {
  const { format, data, indices, indptr, shape } = sparse;
  const [nRows, nCols] = shape;
  const dense = new Float64Array(nRows * nCols);

  if (format === "csr") {
    for (let row = 0; row < nRows; row++) {
      for (let j = indptr[row]; j < indptr[row + 1]; j++) {
        dense[row * nCols + indices[j]] = data[j];
      }
    }
  } else if (format === "csc") {
    for (let col = 0; col < nCols; col++) {
      for (let j = indptr[col]; j < indptr[col + 1]; j++) {
        dense[indices[j] * nCols + col] = data[j];
      }
    }
  } else {
    throw new Error(`Unknown sparse format: "${format}"`);
  }

  return { data: dense, shape };
}

export async function decodeNode(location, open = defaultOpen) {
  const attrs = location.attrs;
  const encodingType = attrs?.["encoding-type"];

  switch (encodingType) {
    case "anndata":
      throw new Error(
        'Use AnnDataStore.open() to read an anndata root, not decodeNode()',
      );
    case "dataframe":
      return decodeDataframe(location, open);
    case "csr_matrix":
    case "csc_matrix":
      return decodeSparseMatrix(location, open);
    case "categorical":
      return decodeCategorical(location, open);
    case "nullable-integer":
    case "nullable-boolean":
      return decodeNullable(location, open);
    default: {
      // Dense array or unknown group — try as array first
      try {
        const arr = await open(location, { kind: "array" });
        return readArray(arr);
      } catch {
        // Return the group attrs for unrecognized groups
        return { attrs };
      }
    }
  }
}
