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
    console.debug("[usePostMessage] Listener registered, allowedOrigin:", allowedOrigin);

    const listener = (event) => {
      if (allowedOrigin !== "*" && event.origin !== allowedOrigin) {
        console.warn("[usePostMessage] Rejected message from origin:", event.origin, "— expected:", allowedOrigin);
        return;
      }

      const { data } = event;
      if (!data || typeof data !== "object" || typeof data.type !== "string") return;

      console.debug("[usePostMessage] Received message:", { type: data.type, origin: event.origin, payload: data.payload });

      const handler = handlers[data.type];
      if (typeof handler === "function") {
        handler(data.payload);
      } else {
        console.warn("[usePostMessage] No handler for message type:", data.type);
      }
    };

    window.addEventListener("message", listener);
    return () => {
      console.debug("[usePostMessage] Listener removed");
      window.removeEventListener("message", listener);
    };
  }, [handlers, allowedOrigin]);
}
