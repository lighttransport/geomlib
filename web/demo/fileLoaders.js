/**
 * File loaders for point clouds and Gaussian splats
 * Supports: PLY (ASCII/Binary), compressed PLY, and SPZ formats
 */

import * as pako from 'pako';

// ============================================================================
// PLY Parser
// ============================================================================

export async function loadPLY(file) {
  const arrayBuffer = await file.arrayBuffer();
  const decoder = new TextDecoder('utf-8');

  // Read header to determine format
  const headerBytes = new Uint8Array(arrayBuffer, 0, Math.min(8192, arrayBuffer.byteLength));
  const headerText = decoder.decode(headerBytes);

  // Check if compressed
  if (headerText.includes('comment compressed')) {
    return parsePLYCompressed(arrayBuffer);
  }

  // Parse standard PLY
  const isAscii = headerText.includes('format ascii');
  const isBinary = headerText.includes('format binary_little_endian') ||
                   headerText.includes('format binary_big_endian');

  if (isAscii) {
    return parsePLYAscii(arrayBuffer);
  } else if (isBinary) {
    return parsePLYBinary(arrayBuffer);
  } else {
    throw new Error('Unknown PLY format');
  }
}

function parsePLYAscii(arrayBuffer) {
  const decoder = new TextDecoder('utf-8');
  const text = decoder.decode(arrayBuffer);
  const lines = text.split('\n');

  let headerEnd = 0;
  let vertexCount = 0;
  let properties = [];
  let inHeader = true;

  // Parse header
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();

    if (line.startsWith('element vertex')) {
      vertexCount = parseInt(line.split(' ')[2]);
    } else if (line.startsWith('property')) {
      const parts = line.split(' ');
      properties.push({ type: parts[1], name: parts[2] });
    } else if (line === 'end_header') {
      headerEnd = i + 1;
      inHeader = false;
      break;
    }
  }

  // Parse vertex data
  const positions = [];
  const normals = [];
  const colors = [];

  const hasNormals = properties.some(p => p.name === 'nx' || p.name === 'normal_x');
  const hasColors = properties.some(p => p.name === 'red' || p.name === 'r');

  for (let i = headerEnd; i < headerEnd + vertexCount && i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line) continue;

    const values = line.split(/\s+/).map(v => parseFloat(v));

    // Assume first 3 are x, y, z
    positions.push(values[0], values[1], values[2]);

    if (hasNormals && values.length >= 6) {
      normals.push(values[3], values[4], values[5]);
    }

    if (hasColors) {
      const colorStart = hasNormals ? 6 : 3;
      if (values.length >= colorStart + 3) {
        colors.push(values[colorStart], values[colorStart + 1], values[colorStart + 2]);
      }
    }
  }

  return {
    positions: new Float32Array(positions),
    normals: normals.length > 0 ? new Float32Array(normals) : null,
    colors: colors.length > 0 ? new Float32Array(colors) : null,
    vertexCount: positions.length / 3
  };
}

function parsePLYBinary(arrayBuffer) {
  const decoder = new TextDecoder('utf-8');
  const headerBytes = new Uint8Array(arrayBuffer, 0, Math.min(8192, arrayBuffer.byteLength));
  const headerText = decoder.decode(headerBytes);
  const lines = headerText.split('\n');

  let headerEnd = 0;
  let vertexCount = 0;
  let properties = [];

  // Parse header
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();

    if (line.startsWith('element vertex')) {
      vertexCount = parseInt(line.split(' ')[2]);
    } else if (line.startsWith('property')) {
      const parts = line.split(' ');
      properties.push({ type: parts[1], name: parts[2] });
    } else if (line === 'end_header') {
      // Find actual byte position of end_header
      const headerStr = lines.slice(0, i + 1).join('\n') + '\n';
      headerEnd = new TextEncoder().encode(headerStr).length;
      break;
    }
  }

  // Parse binary data
  const dataView = new DataView(arrayBuffer, headerEnd);
  const positions = [];
  const normals = [];
  const colors = [];

  const hasNormals = properties.some(p => p.name === 'nx' || p.name === 'normal_x');
  const hasColors = properties.some(p => p.name === 'red' || p.name === 'r');

  let offset = 0;
  for (let i = 0; i < vertexCount; i++) {
    // Read x, y, z (assuming float)
    const x = dataView.getFloat32(offset, true); offset += 4;
    const y = dataView.getFloat32(offset, true); offset += 4;
    const z = dataView.getFloat32(offset, true); offset += 4;
    positions.push(x, y, z);

    // Read normals if present
    if (hasNormals) {
      const nx = dataView.getFloat32(offset, true); offset += 4;
      const ny = dataView.getFloat32(offset, true); offset += 4;
      const nz = dataView.getFloat32(offset, true); offset += 4;
      normals.push(nx, ny, nz);
    }

    // Read colors if present (assuming uchar RGB)
    if (hasColors) {
      const r = dataView.getUint8(offset); offset += 1;
      const g = dataView.getUint8(offset); offset += 1;
      const b = dataView.getUint8(offset); offset += 1;
      colors.push(r / 255, g / 255, b / 255);
    }
  }

  return {
    positions: new Float32Array(positions),
    normals: normals.length > 0 ? new Float32Array(normals) : null,
    colors: colors.length > 0 ? new Float32Array(colors) : null,
    vertexCount: positions.length / 3
  };
}

function parsePLYCompressed(arrayBuffer) {
  // Compressed PLY uses gzip/deflate after header
  const decoder = new TextDecoder('utf-8');
  const headerBytes = new Uint8Array(arrayBuffer, 0, Math.min(8192, arrayBuffer.byteLength));
  const headerText = decoder.decode(headerBytes);
  const lines = headerText.split('\n');

  let headerEnd = 0;
  for (let i = 0; i < lines.length; i++) {
    if (lines[i].trim() === 'end_header') {
      const headerStr = lines.slice(0, i + 1).join('\n') + '\n';
      headerEnd = new TextEncoder().encode(headerStr).length;
      break;
    }
  }

  // Decompress data
  const compressedData = new Uint8Array(arrayBuffer, headerEnd);
  const decompressed = pako.inflate(compressedData);

  // Create new buffer with header + decompressed data
  const headerBuffer = new Uint8Array(arrayBuffer, 0, headerEnd);
  const newBuffer = new Uint8Array(headerEnd + decompressed.length);
  newBuffer.set(headerBuffer);
  newBuffer.set(decompressed, headerEnd);

  // Parse as binary PLY
  return parsePLYBinary(newBuffer.buffer);
}

// ============================================================================
// SPZ Parser (Gaussian Splat Format)
// ============================================================================

export async function loadSPZ(file) {
  const arrayBuffer = await file.arrayBuffer();

  // SPZ is typically compressed
  let data;
  try {
    const compressed = new Uint8Array(arrayBuffer);
    data = pako.inflate(compressed);
  } catch (e) {
    // Maybe not compressed
    data = new Uint8Array(arrayBuffer);
  }

  // Parse Gaussian splat data
  // Format: positions (float32 x3), scales (float32 x3), rotations (float32 x4), colors (float32 x3), opacity (float32)
  const dataView = new DataView(data.buffer);
  const floatSize = 4;
  const splatSize = (3 + 3 + 4 + 3 + 1) * floatSize; // 14 floats per splat
  const splatCount = Math.floor(data.length / splatSize);

  const positions = [];
  const scales = [];
  const rotations = [];
  const colors = [];
  const opacities = [];

  let offset = 0;
  for (let i = 0; i < splatCount; i++) {
    // Position
    positions.push(
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true)
    ); offset += 4;

    // Scale
    scales.push(
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true)
    ); offset += 4;

    // Rotation (quaternion)
    rotations.push(
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true)
    ); offset += 4;

    // Color (SH or RGB)
    colors.push(
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true), offset += 4,
      dataView.getFloat32(offset, true)
    ); offset += 4;

    // Opacity
    opacities.push(dataView.getFloat32(offset, true));
    offset += 4;
  }

  return {
    positions: new Float32Array(positions),
    scales: new Float32Array(scales),
    rotations: new Float32Array(rotations),
    colors: new Float32Array(colors),
    opacities: new Float32Array(opacities),
    splatCount
  };
}

// ============================================================================
// Generic File Loader
// ============================================================================

export async function loadPointCloudFile(file) {
  const ext = file.name.split('.').pop().toLowerCase();

  switch (ext) {
    case 'ply':
      return { type: 'ply', data: await loadPLY(file) };

    case 'spz':
    case 'splat':
      return { type: 'splat', data: await loadSPZ(file) };

    default:
      throw new Error(`Unsupported file format: ${ext}`);
  }
}
