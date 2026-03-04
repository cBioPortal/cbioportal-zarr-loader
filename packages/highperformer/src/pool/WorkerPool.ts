interface PoolWorker {
  worker: Worker
  busy: boolean
}

interface QueuedTask {
  message: Record<string, unknown>
  resolve: (value: unknown) => void
  reject: (reason: unknown) => void
  taskId: number
}

export class WorkerPool {
  private workers: PoolWorker[]
  private queue: QueuedTask[] = []
  private pending = new Map<number, (value: unknown) => void>()
  private nextTaskId = 1

  constructor(factory: () => Worker, size?: number) {
    const poolSize = size ?? Math.max(1, (navigator.hardwareConcurrency ?? 4) - 1)
    this.workers = Array.from({ length: poolSize }, () => {
      const worker = factory()
      const pw: PoolWorker = { worker, busy: false }

      worker.onmessage = (e: MessageEvent) => {
        const { _poolTaskId, ...result } = e.data
        const resolver = this.pending.get(_poolTaskId)
        if (resolver) {
          this.pending.delete(_poolTaskId)
          resolver(result)
        }
        pw.busy = false
        this.dequeue()
      }

      return pw
    })
  }

  dispatch<T = unknown>(message: Record<string, unknown>): Promise<T> {
    const taskId = this.nextTaskId++
    return new Promise<T>((resolve, reject) => {
      const idle = this.workers.find((w) => !w.busy)
      if (idle) {
        this.send(idle, { ...message, _poolTaskId: taskId }, resolve as (v: unknown) => void, taskId)
      } else {
        this.queue.push({ message, resolve: resolve as (v: unknown) => void, reject, taskId })
      }
    })
  }

  dispose(): void {
    for (const pw of this.workers) {
      pw.worker.terminate()
    }
    this.workers = []
    this.queue = []
    this.pending.clear()
  }

  private send(pw: PoolWorker, message: Record<string, unknown>, resolve: (v: unknown) => void, taskId: number): void {
    pw.busy = true
    this.pending.set(taskId, resolve)
    pw.worker.postMessage(message)
  }

  private dequeue(): void {
    if (this.queue.length === 0) return
    const idle = this.workers.find((w) => !w.busy)
    if (!idle) return
    const task = this.queue.shift()!
    this.send(idle, { ...task.message, _poolTaskId: task.taskId }, task.resolve, task.taskId)
  }
}
