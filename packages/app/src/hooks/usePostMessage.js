import { useEffect } from "react";

/**
 * Listen for window postMessage events and dispatch to handlers by type.
 *
 * Messages must follow the envelope format: { type: string, payload: any }
 * Unknown types are silently ignored.
 *
 * @param {Object} handlers - Map of type string to async handler function
 * @param {string} [allowedOrigin="*"] - Origin to accept messages from ("*" = any)
 */
export default function usePostMessage(handlers, allowedOrigin = "*") {
  useEffect(() => {
    const listener = (event) => {
      if (allowedOrigin !== "*" && event.origin !== allowedOrigin) return;

      const { data } = event;
      if (!data || typeof data !== "object" || typeof data.type !== "string") return;

      const handler = handlers[data.type];
      if (typeof handler === "function") {
        handler(data.payload);
      }
    };

    window.addEventListener("message", listener);
    return () => window.removeEventListener("message", listener);
  }, [handlers, allowedOrigin]);
}
