import { z } from "zod";

export const ColorBySchema = z.object({
  type: z.enum(["category", "gene"]),
  value: z.string(),
  color_scale: z.enum(["viridis", "magma"]).optional(),
});

export const SelectionSchema = z.object({
  target: z.string(),
  values: z.array(z.union([z.string(), z.number()])),
});

export const ViewSchema = z.object({
  name: z.string().optional(),
  embedding_key: z.string().optional(),
  selection: SelectionSchema,
  active_tooltips: z.array(z.string()).optional(),
  color_by: ColorBySchema.optional(),
});

export const DefaultsSchema = z.object({
  embedding_key: z.string().optional(),
  active_tooltips: z.array(z.string()).optional(),
  color_by: ColorBySchema.optional(),
});

export const FilterSchema = z.object({
  defaults: DefaultsSchema.optional(),
  initial_view: z.union([z.string(), z.number().int().nonnegative()]),
  saved_views: z.array(ViewSchema),
});
